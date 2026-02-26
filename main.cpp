#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <climits>
#include <fstream>
#include <sstream>

extern "C" {
#include "gwfa.h"
extern int gfa_ed_dbg;
}

/* ---- 2-bit encoding: A=0 C=1 G=2 T=3 ---- */

static inline uint32_t char_to_2bit(char c) {
	switch (c) {
	case 'A': case 'a': return 0;
	case 'C': case 'c': return 1;
	case 'G': case 'g': return 2;
	case 'T': case 't': return 3;
	default: return 0;
	}
}

// Pack a char string into uint32_t array
// (16 chars per word, LSB-first)
static uint32_t *encode_2bit(
	const char *s, size_t len)
{
	size_t nw = (len + 15) / 16;
	uint32_t *out = (uint32_t*)calloc(
		nw, sizeof(uint32_t));
	for (size_t i = 0; i < len; ++i)
		out[i >> 4] |=
			char_to_2bit(s[i]) << ((i & 0xF) << 1);
	return out;
}

/* ---- Dump loading ---- */

struct GwfaIterInput {
	int32_t ql;
	uint32_t *q_enc;  // 2-bit packed query
	int32_t s_term;
	subgfa_subgraph_t *sub;
};

static FILE *openDumpFile(
	const std::string &dir, const char *name)
{
	std::string path = dir + "/" + name;
	FILE *fp = fopen(path.c_str(), "r");
	if (!fp) {
		fprintf(stderr,
			"cannot open %s\n", path.c_str());
		exit(1);
	}
	return fp;
}

static std::vector<uint32_t> parseU32Line(
	const std::string &line)
{
	std::vector<uint32_t> v;
	std::istringstream iss(line);
	uint32_t x;
	while (iss >> x) v.push_back(x);
	return v;
}

static std::vector<int32_t> parseI32Line(
	const std::string &line)
{
	std::vector<int32_t> v;
	std::istringstream iss(line);
	int32_t x;
	while (iss >> x) v.push_back(x);
	return v;
}

static std::vector<uint64_t> parseU64Line(
	const std::string &line)
{
	std::vector<uint64_t> v;
	std::istringstream iss(line);
	uint64_t x;
	while (iss >> x) v.push_back(x);
	return v;
}

static bool readLine(FILE *fp, std::string &out)
{
	char buf[1 << 20];
	if (!fgets(buf, sizeof(buf), fp))
		return false;
	size_t len = strlen(buf);
	while (len > 0
		&& (buf[len-1]=='\n'
		|| buf[len-1]=='\r'))
		--len;
	out.assign(buf, len);
	return true;
}

static subgfa_subgraph_t *buildSubgraph(
	uint32_t n_vtx, uint64_t n_arc,
	const std::string &graphSeq,
	const std::vector<uint32_t> &seq_off,
	const std::vector<int32_t> &seq_len,
	const std::vector<uint32_t> &arc_v,
	const std::vector<uint32_t> &arc_w,
	const std::vector<int32_t> &arc_ow,
	const std::vector<uint64_t> &idx)
{
	subgfa_subgraph_t *sub =
		(subgfa_subgraph_t*)calloc(
			1, sizeof(subgfa_subgraph_t));
	sub->n_vtx = n_vtx;
	sub->n_arc = n_arc;

	sub->graphSeq = encode_2bit(
		graphSeq.data(), graphSeq.size());

	sub->seq_off = (uint32_t*)
		malloc(n_vtx * sizeof(uint32_t));
	memcpy(sub->seq_off, seq_off.data(),
		n_vtx * sizeof(uint32_t));

	sub->seq_len = (int32_t*)
		malloc(n_vtx * sizeof(int32_t));
	memcpy(sub->seq_len, seq_len.data(),
		n_vtx * sizeof(int32_t));

	sub->arc = (subgfa_arc_t*)
		malloc(n_arc * sizeof(subgfa_arc_t));
	for (uint64_t i = 0; i < n_arc; i++) {
		sub->arc[i].v = arc_v[i];
		sub->arc[i].w = arc_w[i];
		sub->arc[i].ow = arc_ow[i];
	}

	sub->idx = (uint64_t*)
		malloc(n_vtx * sizeof(uint64_t));
	memcpy(sub->idx, idx.data(),
		n_vtx * sizeof(uint64_t));

	return sub;
}

static std::vector<GwfaIterInput>
loadGwfaDump(const std::string &dumpDir)
{
	FILE *fql   = openDumpFile(dumpDir, "ql.txt");
	FILE *fq    = openDumpFile(dumpDir, "q.txt");
	FILE *fst   = openDumpFile(
		dumpDir, "s_term.txt");
	FILE *fnv   = openDumpFile(
		dumpDir, "n_vtx.txt");
	FILE *fna   = openDumpFile(
		dumpDir, "n_arc.txt");
	FILE *fgs   = openDumpFile(
		dumpDir, "graphSeq.txt");
	FILE *fsoff = openDumpFile(
		dumpDir, "seq_off.txt");
	FILE *fslen = openDumpFile(
		dumpDir, "seq_len.txt");
	FILE *fav   = openDumpFile(
		dumpDir, "arc_v.txt");
	FILE *faw   = openDumpFile(
		dumpDir, "arc_w.txt");
	FILE *faow  = openDumpFile(
		dumpDir, "arc_ow.txt");
	FILE *fidx  = openDumpFile(
		dumpDir, "idx.txt");

	std::vector<GwfaIterInput> result;
	std::string line;

	while (readLine(fql, line)) {
		GwfaIterInput inp;
		inp.ql = std::stoi(line);

		std::string q_str;
		readLine(fq, q_str);
		inp.q_enc = encode_2bit(
			q_str.data(), q_str.size());

		readLine(fst, line);
		inp.s_term = std::stoi(line);

		/* subgraph */
		readLine(fnv, line);
		uint32_t nv =
			(uint32_t)std::stoul(line);

		if (nv == 0) {
			readLine(fna, line);
			readLine(fgs, line);
			readLine(fsoff, line);
			readLine(fslen, line);
			readLine(fav, line);
			readLine(faw, line);
			readLine(faow, line);
			readLine(fidx, line);
			inp.sub = NULL;
		} else {
			readLine(fna, line);
			uint64_t na = std::stoull(line);

			std::string gs;
			readLine(fgs, gs);

			readLine(fsoff, line);
			auto seq_off = parseU32Line(line);
			readLine(fslen, line);
			auto seq_len = parseI32Line(line);
			readLine(fav, line);
			auto av = parseU32Line(line);
			readLine(faw, line);
			auto aw = parseU32Line(line);
			readLine(faow, line);
			auto aow = parseI32Line(line);
			readLine(fidx, line);
			auto idx = parseU64Line(line);

			inp.sub = buildSubgraph(nv, na,
				gs, seq_off, seq_len,
				av, aw, aow, idx);
		}
		result.push_back(std::move(inp));
	}

	fclose(fql);  fclose(fq);
	fclose(fst);
	fclose(fnv);  fclose(fna);
	fclose(fgs);  fclose(fsoff);
	fclose(fslen); fclose(fav);
	fclose(faw);  fclose(faow);
	fclose(fidx);
	return result;
}

/* ---- Main ---- */

int main(int argc, char *argv[])
{
	std::string dumpDir = "GwfaDump";
	if (argc >= 2) dumpDir = argv[1];
	int numItersLimit = INT_MAX;
	if (argc >= 3)
		numItersLimit = std::stoi(argv[2]);

	std::cout << "Loading Inputs" << std::endl;
	auto load_start =
		std::chrono::system_clock::now();
	auto inputs = loadGwfaDump(dumpDir);
	auto load_end =
		std::chrono::system_clock::now();

	int n = std::min((int)inputs.size(),
		numItersLimit);

	FILE *sfp = fopen("scores.txt", "w");

	std::cout << "Running Kernel" << std::endl;
	auto kernel_start =
		std::chrono::system_clock::now();
	for (int i = 0; i < n; i++) {
		auto &inp = inputs[i];
		auto tS =
			std::chrono::system_clock::now();
		int score;
		if (!inp.sub) {
			score = -1;
		} else {
			score = gwfa(inp.ql,
				inp.q_enc,
				inp.sub, inp.s_term,
				gfa_ed_dbg);
		}
		auto tE =
			std::chrono::system_clock::now();
		auto us = std::chrono::duration_cast<
			std::chrono::microseconds>(
			tE - tS).count();
		fprintf(sfp, "%d\n", score);
		fflush(sfp);
		std::cout << "i: " << i << std::endl;
		std::cout << "iterTime: " << us
			<< "us" << std::endl;
		std::cout << "ql: " << inp.ql
			<< std::endl;
		std::cout << std::endl;
		free(inp.q_enc);
		inp.q_enc = NULL;
		if (inp.sub) {
			subgfa_subgraph_destroy(inp.sub);
			inp.sub = NULL;
		}
	}
	auto kernel_end =
		std::chrono::system_clock::now();
	fclose(sfp);
	std::cout << "Kernel Complete" << std::endl;

	auto load_us =
		std::chrono::duration_cast<
		std::chrono::microseconds>(
		load_end - load_start).count();
	auto kernel_us =
		std::chrono::duration_cast<
		std::chrono::microseconds>(
		kernel_end - kernel_start).count();
	std::cout << "load time: " << load_us
		<< "us" << std::endl;
	std::cout << "kernel time: " << kernel_us
		<< "us" << std::endl;
}
