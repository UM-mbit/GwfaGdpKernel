#ifndef GWFA_H
#define GWFA_H

#include <stdint.h>

typedef struct {
	uint32_t v;
	uint32_t w;
	int32_t ow;
} subgfa_arc_t;

typedef struct {
	uint32_t n_vtx;
	uint64_t n_arc;
	char *graphSeq;
	uint32_t *seq_off;
	int32_t *seq_len;
	subgfa_arc_t *arc;
	uint64_t *idx;
} subgfa_subgraph_t;

#ifdef __cplusplus
extern "C" {
#endif

int gwfa(int32_t ql, const char *q,
	uint32_t startV, uint32_t endV,
	subgfa_subgraph_t *sub, int32_t s_term,
	int dbg);

void subgfa_subgraph_destroy(
	subgfa_subgraph_t *sub);

#ifdef __cplusplus
}
#endif

#endif
