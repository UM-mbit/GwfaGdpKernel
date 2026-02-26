#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include "gwfa.h"
#include "ksort.h"

/* 2-bit packed encoding: A=0 C=1 G=2 T=3,
   16 chars per uint32_t, LSB-first */
#define GET_2BIT(arr, idx) \
	(((arr)[(idx) >> 4] >> (((idx) & 0xF) << 1)) & 0x3)

#ifndef GFA_ED_DBG
#define GFA_ED_DBG 0
#endif
int gfa_ed_dbg = GFA_ED_DBG;

static FILE *gwf_scores_fp(void) {
	static FILE *fp;
	if (!fp) fp = fopen("scores.txt", "a");
	return fp;
}

static FILE *gwf_wf_debug_fp(void) {
	static FILE *fp;
	if (!fp) fp = fopen("wfDebug.txt", "w");
	return fp;
}

/*****************
 * Edit distance *
 *****************/

/* Static array capacity constants */
#define DIAG_CAP  (16 << 20)   /* 16M = 2^24 */
#define INTV_CAP  (1 << 21)    /* 2M (~32 MB) */
#define HA_BITS   22           /* 4M slots */
#define HA_CAP    (1 << HA_BITS)
#define HA_MASK   (HA_CAP - 1)
#define A_MASK    (DIAG_CAP - 1)

#define GWF_DIAG_SHIFT 0x4000

static inline uint32_t gwf_gen_vd(
	uint32_t v, int32_t d)
{
	return (uint32_t)v << 16
		| (GWF_DIAG_SHIFT + d);
}

/*
 * Diagonal interval
 */
typedef struct {
	uint32_t vd0, vd1;
} gwf_intv_t;

#define intvd_key(x) ((x).vd0)
KRADIX_SORT_INIT(gwf_intv, gwf_intv_t,
	intvd_key, 4)

#define subgfa_arc_n(s, v) \
	((uint32_t)(s)->idx[(v)])
#define subgfa_arc_a(s, v) \
	(&(s)->arc[(s)->idx[(v)]>>32])

static int gwf_intv_is_sorted(int32_t n_a,
	const gwf_intv_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].vd0 > a[i].vd0) break;
	return (i == n_a);
}

// merge overlapping intervals; input sorted
static size_t gwf_intv_merge_adj(
	size_t n, gwf_intv_t *a)
{
	size_t i, k;
	uint32_t st, en;
	if (n == 0) return 0;
	st = a[0].vd0, en = a[0].vd1;
	for (i = 1, k = 0; i < n; ++i) {
		if (a[i].vd0 > en) {
			a[k].vd0 = st, a[k++].vd1 = en;
			st = a[i].vd0, en = a[i].vd1;
		} else en = en > a[i].vd1
			? en : a[i].vd1;
	}
	a[k].vd0 = st, a[k++].vd1 = en;
	return k;
}

// merge two sorted interval lists
static size_t gwf_intv_merge2(gwf_intv_t *a,
	size_t n_b, const gwf_intv_t *b,
	size_t n_c, const gwf_intv_t *c)
{
	size_t i = 0, j = 0, k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd0 <= c[j].vd0)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
	return gwf_intv_merge_adj(k, a);
}

/*
 * Diagonal
 */
typedef struct {
	uint32_t vd;
	int32_t k;
} gwf_diag_t;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 4)

/*
 * Core GWFA routine
 */

/* File-scope static arrays */
static gwf_diag_t s_diag_a[DIAG_CAP];
static gwf_diag_t s_diag_b[DIAG_CAP];
static gwf_diag_t s_A[DIAG_CAP];
static gwf_intv_t s_intv[INTV_CAP];
static gwf_intv_t s_tmp[INTV_CAP];
static gwf_intv_t s_swap[INTV_CAP];
static gwf_diag_t s_sort_buf[DIAG_CAP];
static uint32_t   s_ha_keys[HA_CAP];
static uint8_t    s_ha_occ[HA_CAP];
static uint32_t   s_ha_dirty[HA_CAP];

/* File-scope mutable counters */
static uint32_t s_ha_n_dirty;
static uint32_t s_A_head, s_A_tail, s_A_count;
static size_t   s_intv_n, s_tmp_n;

/* ---- hash set helpers ---- */
static inline void ha_clear(void) {
	uint32_t i;
	for (i = 0; i < s_ha_n_dirty; ++i)
		s_ha_occ[s_ha_dirty[i]] = 0;
	s_ha_n_dirty = 0;
}

static inline uint32_t ha_put(
	uint32_t key, int *absent)
{
	uint32_t h = (uint32_t)key * 2654435769U
		>> (32 - HA_BITS);
	while (s_ha_occ[h]) {
		if (s_ha_keys[h] == key) {
			*absent = 0;
			return h;
		}
		h = (h + 1) & HA_MASK;
	}
	s_ha_keys[h] = key;
	s_ha_occ[h] = 1;
	s_ha_dirty[s_ha_n_dirty++] = h;
	*absent = 1;
	return h;
}

/* ---- queue A helpers ---- */
static inline void A_clear(void) {
	s_A_head = s_A_tail = s_A_count = 0;
}

static inline uint32_t A_size(void) {
	return s_A_count;
}

static inline gwf_diag_t *A_pushp(void) {
	gwf_diag_t *p =
		&s_A[s_A_tail++ & A_MASK];
	s_A_count++;
	return p;
}

static inline gwf_diag_t A_shift(void) {
	gwf_diag_t t =
		s_A[s_A_head++ & A_MASK];
	s_A_count--;
	return t;
}

// push (v,d,k) to the end of an array
static inline void gwf_diag_push(
	gwf_diag_t *a, int32_t *n,
	uint32_t v, int32_t d, int32_t k)
{
	gwf_diag_t *p = &a[(*n)++];
	p->vd = gwf_gen_vd(v, d), p->k = k;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(
	gwf_diag_t *p, uint32_t v,
	int32_t d, int32_t k)
{
	uint32_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->k = p->k > k ? p->k : k;
		return 0;
	}
	return 1;
}

// sort using n_sorted as split point
static void gwf_diag_sort(int32_t n_a,
	gwf_diag_t *a, int32_t n_sorted,
	gwf_diag_t *buf)
{
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	n_b = n_sorted;
	n_c = n_a - n_sorted;
	b = buf, c = b + n_b;
	memcpy(b, a, n_b * sizeof(*a));
	memcpy(c, a + n_b, n_c * sizeof(*a));
	radix_sort_gwf_ed(c, c + n_c);

	i = j = k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd <= c[j].vd)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
}

// remove diagonals not on the wavefront
static int32_t gwf_diag_dedup(int32_t n_a,
	gwf_diag_t *a, int32_t n_sorted,
	gwf_diag_t *buf)
{
	int32_t i, n, st;
	if (n_sorted < n_a)
		gwf_diag_sort(n_a, a, n_sorted, buf);
	for (i = 1, st = 0, n = 0; i <= n_a; ++i) {
		if (i == n_a || a[i].vd != a[st].vd) {
			int32_t j, max_j = st;
			if (st + 1 < i)
				for (j = st + 1; j < i; ++j)
					if (a[max_j].k < a[j].k)
						max_j = j;
			a[n++] = a[max_j];
			st = i;
		}
	}
	return n;
}

// use forbidden bands to remove diags
static int32_t gwf_mixed_dedup(int32_t n_a,
	gwf_diag_t *a, int32_t n_b,
	gwf_intv_t *b)
{
	int32_t i = 0, j = 0, k = 0;
	while (i < n_a && j < n_b) {
		if (a[i].vd >= b[j].vd0
			&& a[i].vd < b[j].vd1) ++i;
		else if (a[i].vd >= b[j].vd1) ++j;
		else a[k++] = a[i++];
	}
	while (i < n_a) a[k++] = a[i++];
	return k;
}

// remove diagonals not on the wavefront
static int32_t gwf_dedup(int32_t n_a,
	gwf_diag_t *a, int32_t n_sorted)
{
	if (s_intv_n + s_tmp_n > 0) {
		size_t swap_n;
		if (!gwf_intv_is_sorted(
			s_tmp_n, s_tmp))
			radix_sort_gwf_intv(s_tmp,
				s_tmp + s_tmp_n);
		memcpy(s_swap, s_intv,
			s_intv_n * sizeof(gwf_intv_t));
		swap_n = s_intv_n;
		s_intv_n = gwf_intv_merge2(s_intv,
			swap_n, s_swap, s_tmp_n, s_tmp);
	}
	n_a = gwf_diag_dedup(
		n_a, a, n_sorted, s_sort_buf);
	if (s_intv_n > 0)
		n_a = gwf_mixed_dedup(n_a, a,
			s_intv_n, s_intv);
	return n_a;
}

// reach the wavefront (scalar 2-bit compare)
static inline int32_t gwf_extend1(
	int32_t d, int32_t k,
	int32_t vl, const uint32_t *ts,
	int32_t ts_off,
	int32_t ql, const uint32_t *qs)
{
	int32_t max_k =
		(ql - d < vl ? ql - d : vl) - 1;
	while (k < max_k
		&& GET_2BIT(ts, ts_off + k + 1)
			== GET_2BIT(qs, d + k + 1))
		++k;
	return k;
}

// emit a b-entry, filtering out-of-bounds
static inline void emit_b(
	gwf_diag_t *B_a, int32_t *B_n,
	uint32_t vd, int32_t k,
	int32_t v, int32_t vl, int32_t ql)
{
	int32_t d = (int32_t)(vd & 0xFFFF)
		- GWF_DIAG_SHIFT;
	if (d + k < ql && k < vl) {
		gwf_diag_t *p = &B_a[(*B_n)++];
		p->vd = vd;
		p->k = k;
	} else if (k == vl) {
		if (s_tmp_n >= INTV_CAP) {
			fprintf(stderr,
				"FATAL: s_tmp overflow "
				"(n=%zu, cap=%d)\n",
				s_tmp_n, INTV_CAP);
			exit(1);
		}
		gwf_intv_t *qi = &s_tmp[s_tmp_n++];
		qi->vd0 = gwf_gen_vd(v, d);
		qi->vd1 = qi->vd0 + 1;
	}
}

// check if diagonal hit boundary -> push to A
static inline void boundary_check(
	gwf_diag_t *p, int32_t vl, int32_t ql)
{
	int32_t d = (int32_t)(p->vd & 0xFFFF)
		- GWF_DIAG_SHIFT;
	if (p->k == vl - 1 || d + p->k == ql - 1)
		*A_pushp() = *p;
}

// wfa_extend and wfa_next combined (fused)
static gwf_diag_t *gwf_ed_extend(
	const subgfa_subgraph_t *sub,
	int32_t s, int32_t ql, const uint32_t *q,
	int32_t *n_a_, gwf_diag_t *a,
	int *terminate)
{
	int32_t i, n = *n_a_, do_dedup = 1;
	gwf_diag_t *B_a, *b;
	int32_t B_n;

	/* B is the OTHER static buffer (ping-pong) */
	B_a = (a == s_diag_a) ? s_diag_b : s_diag_a;
	B_n = 0;

	/* Fused extend + next: process one diagonal
	   at a time, emitting b[] entries inline */
	i = 0;
	while (i < n) {
		int32_t v = a[i].vd >> 16;
		int32_t vl = sub->seq_len[v];
		int32_t ts_off = sub->seq_off[v];
		int32_t ppk, pk, k;

		/* Extend first diagonal, emit b[0] */
		a[i].k = gwf_extend1(
			(int32_t)(a[i].vd & 0xFFFF)
				- GWF_DIAG_SHIFT,
			a[i].k, vl, sub->graphSeq,
			ts_off, ql, q);
		emit_b(B_a, &B_n, a[i].vd - 1,
			a[i].k + 1, v, vl, ql);
		boundary_check(&a[i], vl, ql);
		ppk = a[i].k;
		i++;

		if (i < n
			&& a[i].vd == a[i-1].vd + 1) {
			/* Extend second diagonal,
			   emit b[1] */
			a[i].k = gwf_extend1(
				(int32_t)(a[i].vd & 0xFFFF)
					- GWF_DIAG_SHIFT,
				a[i].k, vl, sub->graphSeq,
				ts_off, ql, q);
			k = (ppk > a[i].k
				? ppk : a[i].k) + 1;
			emit_b(B_a, &B_n, a[i-1].vd,
				k, v, vl, ql);
			boundary_check(&a[i], vl, ql);
			pk = ppk;
			ppk = a[i].k;
			i++;

			/* Interior diagonals */
			while (i < n
				&& a[i].vd == a[i-1].vd + 1)
			{
				a[i].k = gwf_extend1(
					(int32_t)(a[i].vd & 0xFFFF)
						- GWF_DIAG_SHIFT,
					a[i].k, vl,
					sub->graphSeq,
					ts_off, ql, q);
				k = pk;
				if (ppk + 1 > k) k = ppk + 1;
				if (a[i].k + 1 > k)
					k = a[i].k + 1;
				emit_b(B_a, &B_n,
					a[i-1].vd, k,
					v, vl, ql);
				boundary_check(&a[i], vl, ql);
				pk = ppk;
				ppk = a[i].k;
				i++;
			}

			/* Right side: b[n] */
			k = pk > ppk + 1 ? pk : ppk + 1;
			emit_b(B_a, &B_n, a[i-1].vd,
				k, v, vl, ql);
		} else {
			/* n==1: emit b[1] */
			emit_b(B_a, &B_n, a[i-1].vd,
				ppk + 1, v, vl, ql);
		}

		/* Right boundary: b[n+1] */
		emit_b(B_a, &B_n, a[i-1].vd + 1,
			a[i-1].k, v, vl, ql);
	}
	if (A_size() == 0) do_dedup = 0;
	int32_t n_sorted = B_n;

	while (A_size()) {
		gwf_diag_t t;
		uint32_t v;
		int32_t d, k, i, vl;

		t = A_shift();
		v = t.vd >> 16;
		d = (int32_t)(t.vd & 0xFFFF) - GWF_DIAG_SHIFT;
		k = t.k;
		vl = sub->seq_len[v];
		k = gwf_extend1(d, k, vl,
			sub->graphSeq, sub->seq_off[v],
			ql, q);
		i = k + d;

		if (k + 1 < vl && i + 1 < ql) {
			gwf_diag_push(B_a, &B_n, v, d-1, k+1);
			gwf_diag_push(B_a, &B_n, v, d, k+1);
			gwf_diag_push(B_a, &B_n, v, d+1, k);
		} else if (i + 1 < ql) {
			int32_t nv =
				subgfa_arc_n(sub, v);
			int32_t j, n_ext = 0;
			const subgfa_arc_t *av =
				subgfa_arc_a(sub, v);
			gwf_intv_t *p;
			if (s_tmp_n >= INTV_CAP) {
				fprintf(stderr,
					"FATAL: s_tmp overflow in "
					"gwf_ed_extend "
					"(n=%zu, cap=%d)\n",
					s_tmp_n, INTV_CAP);
				exit(1);
			}
			p = &s_tmp[s_tmp_n++];
			p->vd0 = gwf_gen_vd(v, d);
			p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) {
				uint32_t w = av[j].w;
				int32_t ol = av[j].ow;
				int absent;
				ha_put(
					(uint32_t)w<<16
					| ((i + 1) & 0xFFFF),
					&absent);
				if (GET_2BIT(q, i + 1)
				== GET_2BIT(sub->graphSeq,
					sub->seq_off[w] + ol))
				{
					++n_ext;
					if (absent) {
						gwf_diag_t *dp;
						dp = A_pushp();
						dp->vd = gwf_gen_vd(
							w, i + 1 - ol);
						dp->k = ol;
					}
				} else if (absent) {
					gwf_diag_push(B_a, &B_n,
						w, i - ol, ol);
					gwf_diag_push(B_a, &B_n,
						w, i + 1 - ol, ol);
				}
			}
			if (nv == 0 || n_ext != nv)
				gwf_diag_push(
					B_a, &B_n, v, d+1, k);
		} else if (v == GWFA_END_V && k + 1 == vl) {
			*terminate = 1;
			return 0;
		} else if (k + 1 < vl) {
			gwf_diag_push(B_a, &B_n, v, d-1, k+1);
		} else {
			int32_t nv =
				subgfa_arc_n(sub, v), j;
			const subgfa_arc_t *av =
				subgfa_arc_a(sub, v);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(B_a, &B_n,
					av[j].w,
					i - av[j].ow,
					av[j].ow);
		}
	}

	*n_a_ = n = B_n, b = B_a;

	if (do_dedup)
		*n_a_ = n = gwf_dedup(
			n, b, n_sorted);
	return b;
}

static void gwf_ed_print_intv(size_t n,
	gwf_intv_t *a)
{
	FILE *fp = gwf_wf_debug_fp();
	size_t i;
	for (i = 0; i < n; ++i)
		fprintf(fp, "Z\t%d\t%d\t%d\n",
			(int32_t)(a[i].vd0>>16),
			(int32_t)(a[i].vd0 & 0xFFFF)
				- GWF_DIAG_SHIFT,
			(int32_t)(a[i].vd1 & 0xFFFF)
				- GWF_DIAG_SHIFT);
	fflush(fp);
}

static void gwf_ed_print_wf(int32_t n,
	const gwf_diag_t *a)
{
	FILE *fp = gwf_wf_debug_fp();
	int32_t i;
	for (i = 0; i < n; ++i) {
		int32_t nid =
			(int32_t)(a[i].vd >> 16);
		int32_t diag =
			(int32_t)(a[i].vd & 0xFFFF)
			- GWF_DIAG_SHIFT;
		fprintf(fp, "WF\t%d\t%d\t%d\n",
			nid, diag, a[i].k);
	}
	fflush(fp);
}

/* ---- Split API: persistent state ---- */
static subgfa_subgraph_t s_sub_copy;
static int32_t s_ql;
static const uint32_t *s_q;
static int s_dbg;
static gwf_diag_t *s_a;
static int32_t s_n_a;
static int32_t s_last_score;

void gwfa_init(int32_t ql, const uint32_t *q,
	const subgfa_subgraph_t *sub, int dbg)
{
	s_sub_copy = *sub;  /* shallow copy */
	s_ql = ql;
	s_q = q;
	s_dbg = dbg;
	s_last_score = -1;

	s_intv_n = 0;
	s_tmp_n = 0;
	memset(s_ha_occ, 0, sizeof(s_ha_occ));
	s_ha_n_dirty = 0;

	s_a = s_diag_a;
	s_n_a = 1;
	s_a[0].vd = gwf_gen_vd(GWFA_START_V, 0);
	s_a[0].k = -1;
}

void gwfa_reset_step(void)
{
	s_tmp_n = 0;
	ha_clear();
	A_clear();
}

int gwfa_extend_step(int32_t s)
{
	int terminate = 0;
	s_a = gwf_ed_extend(
		&s_sub_copy, s, s_ql, s_q,
		&s_n_a, s_a, &terminate);
	if (terminate) {
		s_last_score = s;
		return 1;
	}
	if (s_n_a == 0) {
		fprintf(stderr,
			"gwfa_extend_step: n_a==0 at s=%d\n",
			s);
		s_last_score = -1;
		return 1;
	}
	return 0;
}

void gwfa_debug_step(int32_t s)
{
	if (s_dbg >= 1) {
		FILE *fp = gwf_wf_debug_fp();
		fprintf(fp,
			"[gfa_ed_step] dist=%d, n=%d,"
			" n_intv=%zd, n_tb=0\n",
			s, s_n_a, s_intv_n);
		fflush(fp);
		gwf_ed_print_intv(s_intv_n, s_intv);
		gwf_ed_print_wf(s_n_a, s_a);
	}
}

int gwfa_get_score(void)
{
	return s_last_score;
}

/* ---- Original monolithic API ---- */

int gwfa(int32_t ql, const uint32_t *q,
	subgfa_subgraph_t *sub, int32_t s_term,
	int dbg)
{
	int32_t s;
	gwfa_init(ql, q, sub, dbg);
	for (s = 0; ; s++) {
		gwfa_reset_step();
		if (gwfa_extend_step(s))
			break;
		if (s >= s_term)
			break;
		gwfa_debug_step(s + 1);
	}
	return gwfa_get_score();
}

void subgfa_subgraph_destroy(
	subgfa_subgraph_t *sub)
{
	if (!sub) return;
	free(sub->graphSeq);
	free(sub->seq_off);
	free(sub->seq_len);
	free(sub->arc);
	free(sub->idx);
	free(sub);
}
