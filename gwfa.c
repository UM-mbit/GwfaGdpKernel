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

/* GWF_DIAG_SHIFT now in gwfa.h */

static inline uint32_t gwf_gen_vd(uint32_t v, int32_t d)
{
	return (uint32_t)v << 16 | (GWF_DIAG_SHIFT + d);
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
	((s)->arc_off[(v)+1] - (s)->arc_off[(v)])
#define subgfa_arc_a(s, v) \
	(&(s)->arc[(s)->arc_off[(v)]])

static int gwf_intv_is_sorted(int32_t n_a, const gwf_intv_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].vd0 > a[i].vd0) break;
	return (i == n_a);
}

// merge overlapping intervals; input sorted
static size_t gwf_intv_merge_adj(size_t n, gwf_intv_t *a)
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
static size_t gwf_intv_merge2(gwf_intv_t *a, size_t n_b, const gwf_intv_t *b,
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

/* MM (main memory) block: all large arrays in one contiguous allocation.
 * Word offsets into s_mm (each entry = 1 int = 4 bytes). */
#define MM_DIAG_A_OFF     0
#define MM_DIAG_B_OFF     (MM_DIAG_A_OFF + DIAG_CAP * 2)
#define MM_A_OFF          (MM_DIAG_B_OFF + DIAG_CAP * 2)
#define MM_INTV_OFF       (MM_A_OFF + DIAG_CAP * 2)
#define MM_NEXT_INTV_OFF  (MM_INTV_OFF + INTV_CAP * 2)
#define MM_SWAP_OFF       (MM_NEXT_INTV_OFF + INTV_CAP * 2)
#define MM_SORT_BUF_OFF   (MM_SWAP_OFF + INTV_CAP * 2)

/* Bucket hash table in MM. Each bucket = 4 keys + 1 count word.
 * Linear probing at bucket granularity. */
#define HA_BUCKET_SIZE  4
#define HA_BUCKET_CAP   (HA_CAP / HA_BUCKET_SIZE)  /* 1M buckets */
#define HA_BUCKET_MASK  (HA_BUCKET_CAP - 1)
#define HA_BUCKET_WORDS 4  /* 4 keys, sentinel-based (no count) */
#define HA_SENTINEL     ((int)0xFFFFFFFF)
#define MM_HA_OFF       (MM_SORT_BUF_OFF + DIAG_CAP * 2)
#define MM_HA_DIRTY_OFF (MM_HA_OFF + HA_BUCKET_CAP * HA_BUCKET_WORDS)
#define MM_TOTAL_WORDS  (MM_HA_DIRTY_OFF + HA_BUCKET_CAP)

static int *s_mm;

static void mm_init(void)
{
	if (s_mm) return;
	s_mm = (int *)malloc((size_t)MM_TOTAL_WORDS * sizeof(int));
	if (!s_mm) {
		fprintf(stderr, "FATAL: cannot allocate MM block (%zu MB)\n",
			(size_t)MM_TOTAL_WORDS * 4 / (1024*1024));
		exit(1);
	}
	/* Fill all bucket slots with sentinel (empty marker) */
	for (size_t i = 0; i < HA_BUCKET_CAP * HA_BUCKET_WORDS; i++)
		s_mm[MM_HA_OFF + i] = HA_SENTINEL;
}

/* Pointer aliases into MM — transparent to existing code */
#define s_diag_a  ((gwf_diag_t*)(s_mm + MM_DIAG_A_OFF))
#define s_diag_b  ((gwf_diag_t*)(s_mm + MM_DIAG_B_OFF))
#define s_A       ((gwf_diag_t*)(s_mm + MM_A_OFF))
#define s_intv    ((gwf_intv_t*)(s_mm + MM_INTV_OFF))
#define s_next_intv_buf \
	((gwf_intv_t*)(s_mm + MM_NEXT_INTV_OFF))
#define s_swap    ((gwf_intv_t*)(s_mm + MM_SWAP_OFF))
#define s_sort_buf ((gwf_diag_t*)(s_mm + MM_SORT_BUF_OFF))

/* File-scope mutable counters */
static uint32_t s_ha_n_dirty;
static uint32_t s_A_head, s_A_tail, s_A_count;
static size_t   s_intv_n, s_next_intv_buf_n;

/* ---- bucket hash set helpers ---- */
/* Clear dirty buckets by writing sentinels to all 4 slots */
static inline void ha_clear(void) {
	for (uint32_t i = 0; i < s_ha_n_dirty; ++i) {
		int b = s_mm[MM_HA_DIRTY_OFF + i];
		int *bp = s_mm + MM_HA_OFF + b * HA_BUCKET_WORDS;
		bp[0] = bp[1] = bp[2] = bp[3] = HA_SENTINEL;
	}
	s_ha_n_dirty = 0;
}

/* Insert key into bucket hash table. Sets *absent=1 if new.
 * Sentinel-based: scan for key match, then sentinel (empty).
 * Linear probing at bucket granularity. */
static inline uint32_t ha_put(uint32_t key, int *absent)
{
	uint32_t h = (uint32_t)key * 2654435769U
		>> (32 - HA_BITS);
	uint32_t b = (h >> 2) & HA_BUCKET_MASK;
	while (1) {
		int *bp = s_mm + MM_HA_OFF + b * HA_BUCKET_WORDS;
		for (int i = 0; i < HA_BUCKET_SIZE; i++) {
			if ((uint32_t)bp[i] == key) {
				*absent = 0;
				return b;
			}
			if (bp[i] == HA_SENTINEL) {
				bp[i] = (int)key;
				if (i == 0) {
					if (s_ha_n_dirty >= HA_BUCKET_CAP) {
						fprintf(stderr, "FATAL: ha_dirty overflow "
							"(n=%u, cap=%d)\n",
							s_ha_n_dirty, HA_BUCKET_CAP);
						exit(1);
					}
					s_mm[MM_HA_DIRTY_OFF + s_ha_n_dirty++]
						= (int)b;
				}
				*absent = 1;
				return b;
			}
		}
		b = (b + 1) & HA_BUCKET_MASK;
	}
}

/* Public wrapper for ha_put (called from pe_array.cpp) */
uint32_t gwfa_ha_put(uint32_t key, int *absent)
{
	return ha_put(key, absent);
}

/* ---- queue A helpers ---- */
static inline void A_clear(void) {
	s_A_head = s_A_tail = s_A_count = 0;
}

static inline uint32_t A_size(void) {
	return s_A_count;
}

static inline gwf_diag_t *A_pushp(void) {
	gwf_diag_t *p = &s_A[s_A_tail++ & A_MASK];
	s_A_count++;
	return p;
}

static inline gwf_diag_t A_shift(void) {
	gwf_diag_t t = s_A[s_A_head++ & A_MASK];
	s_A_count--;
	return t;
}

// push (v,d,k) to the end of an array
static inline void gwf_diag_push(gwf_diag_t *a, int32_t *n, uint32_t v, int32_t d, int32_t k)
{
	gwf_diag_t *p = &a[(*n)++];
	p->vd = gwf_gen_vd(v, d), p->k = k;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(gwf_diag_t *p, uint32_t v, int32_t d, int32_t k)
{
	uint32_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->k = p->k > k ? p->k : k;
		return 0;
	}
	return 1;
}

// sort using n_sorted as split point
static void gwf_diag_sort(int32_t n_a, gwf_diag_t *a, int32_t n_sorted, gwf_diag_t *buf)
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
static int32_t gwf_diag_dedup(int32_t n_a, gwf_diag_t *a, int32_t n_sorted, gwf_diag_t *buf)
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
static int32_t gwf_mixed_dedup(int32_t n_a, gwf_diag_t *a, int32_t n_b, gwf_intv_t *b)
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
static int32_t gwf_dedup(int32_t n_a, gwf_diag_t *a, int32_t n_sorted)
{
	if (s_intv_n + s_next_intv_buf_n > 0) {
		size_t swap_n;
		if (!gwf_intv_is_sorted(s_next_intv_buf_n, s_next_intv_buf))
			radix_sort_gwf_intv(s_next_intv_buf, s_next_intv_buf + s_next_intv_buf_n);
		memcpy(s_swap, s_intv, s_intv_n * sizeof(gwf_intv_t));
		swap_n = s_intv_n;
		s_intv_n = gwf_intv_merge2(s_intv, swap_n, s_swap, s_next_intv_buf_n, s_next_intv_buf);
	}
	n_a = gwf_diag_dedup(n_a, a, n_sorted, s_sort_buf);
	if (s_intv_n > 0)
		n_a = gwf_mixed_dedup(n_a, a, s_intv_n, s_intv);
	return n_a;
}

// reach the wavefront (scalar 2-bit compare)
static inline int32_t gwf_extend1(int32_t d, int32_t k, int32_t vl,
	const uint32_t *ts, int32_t ts_off, int32_t ql, const uint32_t *qs)
{
	int32_t max_k = (ql - d < vl ? ql - d : vl) - 1;
	while (k < max_k && GET_2BIT(ts, ts_off + k + 1) == GET_2BIT(qs, d + k + 1))
		++k;
	return k;
}

// emit a b-entry, filtering out-of-bounds
static inline void emit_b(gwf_diag_t *B_a, int32_t *B_n, uint32_t vd, int32_t k,
	int32_t v, int32_t vl, int32_t ql)
{
	int32_t d = (int32_t)(vd & 0xFFFF) - GWF_DIAG_SHIFT;
	if (d + k < ql && k < vl) {
		gwf_diag_t *p = &B_a[(*B_n)++];
		p->vd = vd;
		p->k = k;
	} else if (k == vl) {
		if (s_next_intv_buf_n >= INTV_CAP) {
			fprintf(stderr, "FATAL: s_next_intv_buf overflow (n=%zu, cap=%d)\n",
				s_next_intv_buf_n, INTV_CAP);
			exit(1);
		}
		gwf_intv_t *qi = &s_next_intv_buf[s_next_intv_buf_n++];
		qi->vd0 = gwf_gen_vd(v, d);
		qi->vd1 = qi->vd0 + 1;
	}
}

// emit_b variant that writes intervals to SPM tile instead of globals
static inline void emit_b_tile(gwf_diag_t *B_a, int32_t *B_n,
	int *tile_intv, int32_t *tile_intv_n, uint32_t vd, int32_t k,
	int32_t v, int32_t vl, int32_t ql)
{
	int32_t d = (int32_t)(vd & 0xFFFF) - GWF_DIAG_SHIFT;
	if (d + k < ql && k < vl) {
		gwf_diag_t *p = &B_a[(*B_n)++];
		p->vd = vd;
		p->k = k;
	} else if (k == vl) {
		uint32_t vd0 = gwf_gen_vd(v, d);
		tile_intv[2 * (*tile_intv_n)] = (int)vd0;
		tile_intv[2 * (*tile_intv_n) + 1] = (int)(vd0 + 1);
		(*tile_intv_n)++;
	}
}

// check if diagonal hit boundary -> push to A
static inline void boundary_check(gwf_diag_t *p, int32_t vl, int32_t ql)
{
	int32_t d = (int32_t)(p->vd & 0xFFFF) - GWF_DIAG_SHIFT;
	if (p->k == vl - 1 || d + p->k == ql - 1)
		*A_pushp() = *p;
}

// wfa_extend and wfa_next combined (fused)
static gwf_diag_t *gwf_ed_extend(const subgfa_subgraph_t *sub,
	int32_t s, int32_t ql, const uint32_t *q,
	int32_t *n_a_, gwf_diag_t *a, int *terminate)
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
		a[i].k = gwf_extend1((int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
			a[i].k, vl, sub->graphSeq, ts_off, ql, q);
		emit_b(B_a, &B_n, a[i].vd - 1, a[i].k + 1, v, vl, ql);
		boundary_check(&a[i], vl, ql);
		ppk = a[i].k;
		i++;

		if (i < n && a[i].vd == a[i-1].vd + 1) {
			/* Extend second diagonal, emit b[1] */
			a[i].k = gwf_extend1((int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
				a[i].k, vl, sub->graphSeq, ts_off, ql, q);
			k = (ppk > a[i].k ? ppk : a[i].k) + 1;
			emit_b(B_a, &B_n, a[i-1].vd, k, v, vl, ql);
			boundary_check(&a[i], vl, ql);
			pk = ppk;
			ppk = a[i].k;
			i++;

			/* Interior diagonals */
			while (i < n && a[i].vd == a[i-1].vd + 1) {
				a[i].k = gwf_extend1((int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
					a[i].k, vl, sub->graphSeq, ts_off, ql, q);
				k = pk;
				if (ppk + 1 > k) k = ppk + 1;
				if (a[i].k + 1 > k) k = a[i].k + 1;
				emit_b(B_a, &B_n, a[i-1].vd, k, v, vl, ql);
				boundary_check(&a[i], vl, ql);
				pk = ppk;
				ppk = a[i].k;
				i++;
			}

			/* Right side: b[n] */
			k = pk > ppk + 1 ? pk : ppk + 1;
			emit_b(B_a, &B_n, a[i-1].vd, k, v, vl, ql);
		} else {
			/* n==1: emit b[1] */
			emit_b(B_a, &B_n, a[i-1].vd, ppk + 1, v, vl, ql);
		}

		/* Right boundary: b[n+1] */
		emit_b(B_a, &B_n, a[i-1].vd + 1, a[i-1].k, v, vl, ql);
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
		k = gwf_extend1(d, k, vl, sub->graphSeq, sub->seq_off[v], ql, q);
		i = k + d;

		if (k + 1 < vl && i + 1 < ql) {
			gwf_diag_push(B_a, &B_n, v, d-1, k+1);
			gwf_diag_push(B_a, &B_n, v, d, k+1);
			gwf_diag_push(B_a, &B_n, v, d+1, k);
		} else if (i + 1 < ql) {
			int32_t nv = subgfa_arc_n(sub, v);
			int32_t j, n_ext = 0;
			const subgfa_arc_t *av = subgfa_arc_a(sub, v);
			gwf_intv_t *p;
			if (s_next_intv_buf_n >= INTV_CAP) {
				fprintf(stderr, "FATAL: s_next_intv_buf overflow in gwf_ed_extend (n=%zu, cap=%d)\n",
					s_next_intv_buf_n, INTV_CAP);
				exit(1);
			}
			p = &s_next_intv_buf[s_next_intv_buf_n++];
			p->vd0 = gwf_gen_vd(v, d);
			p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) {
				uint32_t w = av[j].w;
				int32_t ol = av[j].ow;
				int absent;
				ha_put((uint32_t)w<<16 | ((i + 1) & 0xFFFF), &absent);
				if (GET_2BIT(q, i + 1) == GET_2BIT(sub->graphSeq, sub->seq_off[w] + ol)) {
					++n_ext;
					if (absent) {
						gwf_diag_t *dp = A_pushp();
						dp->vd = gwf_gen_vd(w, i + 1 - ol);
						dp->k = ol;
					}
				} else if (absent) {
					gwf_diag_push(B_a, &B_n, w, i - ol, ol);
					gwf_diag_push(B_a, &B_n, w, i + 1 - ol, ol);
				}
			}
			if (nv == 0 || n_ext != nv)
				gwf_diag_push(B_a, &B_n, v, d+1, k);
		} else if (v == GWFA_END_V && k + 1 == vl) {
			*terminate = 1;
			return 0;
		} else if (k + 1 < vl) {
			gwf_diag_push(B_a, &B_n, v, d-1, k+1);
		} else {
			int32_t nv = subgfa_arc_n(sub, v), j;
			const subgfa_arc_t *av = subgfa_arc_a(sub, v);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(B_a, &B_n, av[j].w, i - av[j].ow, av[j].ow);
		}
	}

	*n_a_ = n = B_n, b = B_a;

	if (do_dedup)
		*n_a_ = n = gwf_dedup(n, b, n_sorted);
	return b;
}

static void gwf_ed_print_intv(size_t n, gwf_intv_t *a)
{
	FILE *fp = gwf_wf_debug_fp();
	size_t i;
	for (i = 0; i < n; ++i)
		fprintf(fp, "Z\t%d\t%d\t%d\n", (int32_t)(a[i].vd0>>16),
			(int32_t)(a[i].vd0 & 0xFFFF) - GWF_DIAG_SHIFT,
			(int32_t)(a[i].vd1 & 0xFFFF) - GWF_DIAG_SHIFT);
	fflush(fp);
}

static void gwf_ed_print_wf(int32_t n, const gwf_diag_t *a)
{
	FILE *fp = gwf_wf_debug_fp();
	int32_t i;
	for (i = 0; i < n; ++i) {
		int32_t nid = (int32_t)(a[i].vd >> 16);
		int32_t diag = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
		fprintf(fp, "WF\t%d\t%d\t%d\n", nid, diag, a[i].k);
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
	mm_init();
	s_sub_copy = *sub;  /* shallow copy */
	s_ql = ql;
	s_q = q;
	s_dbg = dbg;
	s_last_score = -1;

	s_intv_n = 0;
	s_next_intv_buf_n = 0;
	ha_clear();

	A_clear();
	s_a = s_diag_a;
	s_n_a = 1;
	s_a[0].vd = gwf_gen_vd(GWFA_START_V, 0);
	s_a[0].k = -1;
}

void gwfa_reset_step(void)
{
	s_next_intv_buf_n = 0;
	ha_clear();
	A_clear();
}

int gwfa_extend_step(int32_t s)
{
	int terminate = 0;
	s_a = gwf_ed_extend(&s_sub_copy, s, s_ql, s_q, &s_n_a, s_a, &terminate);
	if (terminate) {
		s_last_score = s;
		return 1;
	}
	if (s_n_a == 0) {
		fprintf(stderr, "gwfa_extend_step: n_a==0 at s=%d\n", s);
		s_last_score = -1;
		return 1;
	}
	return 0;
}

/* SPM tile layout constants */
#define N_TILE_DIAGS  64
#define A_TILE_OFF    0     /* 128 words: 64 (vd,k) */
#define SEQ_INFO_OFF  128   /* 128 words: 64 (off,len) */
#define B_TILE_OFF    256   /* 384 words: 192 (vd,k) */
#define INTV_TILE_OFF 640   /* 384 words: 192 intervals (vd0,vd1) */
#define A_OUT_OFF     1024  /* 128 words: 64 (vd,k) */
#define META_OFF      1152  /* 16 words: n_b,n_A,n_seq,tile_n,n_vtx,gs_nw,ql,n_intv */
#define SEQ_REGION_OFF 1280 /* sequences: seq_off, seq_len, graphSeq, query */
#define PE_SPM_SIZE    8192 /* must match PE_SPM_SIZE in sys_def.h */

/* Phase 2 tile SPM layout constants */
#define P2_TILE_SIZE   64
#define P2_INPUT_OFF   0     /* 256 words: 64 × (vd, k, ts_off, vl) */
#define P2_PUSHED_OFF  256   /* 384 words: up to 192 × (vd, k) */
#define P2_INTV_OFF    640   /* 128 words: up to 64 × (vd0, vd1) */
#define P2_FIN0_OFF    768   /* 128 words: up to 64 × (vd, k) */
#define P2_FIN1_OFF    896   /* 128 words: up to 64 × (vd, k) */
#define P2_META_OFF    1024  /* 8 words: metadata */
#define P2_M_PUSHED    0
#define P2_M_INTV      1
#define P2_M_FIN0      2
#define P2_M_FIN1      3
#define P2_M_TILE_N    4

/* Fused extend+next on a slice of diags.
   Reads a[0:n], writes to b[*bn:] and aout[*aoutn:].
   All diags within a group (same vertex, consecutive vd)
   must be in the same call — never split groups.
   seq_off/seq_len/graphSeq/q/ql are explicit so callers can
   pass either main-memory or SPM-resident pointers. */
static void phase1_fused(
	gwf_diag_t *a, int32_t n,
	gwf_diag_t *b, int32_t *bn,
	gwf_diag_t *aout, int32_t *aoutn,
	const uint32_t *seq_off, const int32_t *seq_len,
	const uint32_t *graphSeq, const uint32_t *q, int32_t ql)
{
	int32_t i = 0;
	while (i < n) {
		int32_t v = a[i].vd >> 16;
		int32_t vl = seq_len[v];
		int32_t ts_off = seq_off[v];
		int32_t ppk, pk, k, d;

		a[i].k = gwf_extend1(
			(int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
			a[i].k, vl, graphSeq, ts_off, ql, q);
		emit_b(b, bn, a[i].vd - 1, a[i].k + 1, v, vl, ql);
		d = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
		if (a[i].k == vl - 1 || d + a[i].k == ql - 1)
			aout[(*aoutn)++] = a[i];
		ppk = a[i].k;
		i++;

		if (i < n && a[i].vd == a[i-1].vd + 1) {
			a[i].k = gwf_extend1(
				(int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
				a[i].k, vl, graphSeq, ts_off, ql, q);
			k = (ppk > a[i].k ? ppk : a[i].k) + 1;
			emit_b(b, bn, a[i-1].vd, k, v, vl, ql);
			d = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
			if (a[i].k == vl - 1 || d + a[i].k == ql - 1)
				aout[(*aoutn)++] = a[i];
			pk = ppk;
			ppk = a[i].k;
			i++;

			while (i < n && a[i].vd == a[i-1].vd + 1) {
				a[i].k = gwf_extend1(
					(int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
					a[i].k, vl, graphSeq, ts_off, ql, q);
				k = pk;
				if (ppk + 1 > k) k = ppk + 1;
				if (a[i].k + 1 > k) k = a[i].k + 1;
				emit_b(b, bn, a[i-1].vd, k, v, vl, ql);
				d = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
				if (a[i].k == vl - 1 || d + a[i].k == ql - 1)
					aout[(*aoutn)++] = a[i];
				pk = ppk;
				ppk = a[i].k;
				i++;
			}

			k = pk > ppk + 1 ? pk : ppk + 1;
			emit_b(b, bn, a[i-1].vd, k, v, vl, ql);
		} else {
			emit_b(b, bn, a[i-1].vd, ppk + 1, v, vl, ql);
		}

		emit_b(b, bn, a[i-1].vd + 1, a[i-1].k, v, vl, ql);
	}
}


int gwfa_extend_step_tiled(int32_t s, int *spm)
{
	int32_t cursor;
	gwf_diag_t *B_a;
	int32_t B_n = 0;

	B_a = (s_a == s_diag_a) ? s_diag_b : s_diag_a;

	/* Phase 1: tiled fused extend+next (fixed 64-diag tiles) */
	for (cursor = 0; cursor < s_n_a; ) {
		int32_t remaining = s_n_a - cursor;
		int32_t tile_n = remaining < N_TILE_DIAGS ? remaining : N_TILE_DIAGS;

		/* -- tile_load -- */
		gwf_diag_t *tile_a = (gwf_diag_t *)(spm + A_TILE_OFF);
		memcpy(tile_a, &s_a[cursor], tile_n * sizeof(gwf_diag_t));

		/* -- tile_compute -- */
		gwf_diag_t *tile_b = (gwf_diag_t *)(spm + B_TILE_OFF);
		gwf_diag_t *tile_aout = (gwf_diag_t *)(spm + A_OUT_OFF);
		int32_t tb_n = 0, ta_n = 0;
		phase1_fused(tile_a, tile_n, tile_b, &tb_n, tile_aout, &ta_n,
			s_sub_copy.seq_off, s_sub_copy.seq_len,
			s_sub_copy.graphSeq, s_q, s_ql);

		spm[META_OFF] = tb_n;
		spm[META_OFF + 1] = ta_n;

		/* -- tile_writeback -- */
		memcpy(&B_a[B_n], tile_b, tb_n * sizeof(gwf_diag_t));
		B_n += tb_n;
		for (int32_t j = 0; j < ta_n; j++)
			*A_pushp() = tile_aout[j];

		cursor += tile_n;
	}

	/* Phase 2: A queue processing (same as gwf_ed_extend) */
	int do_dedup = (A_size() > 0);
	int32_t n_sorted = 0; /* tile boundaries may break sort order */

	while (A_size()) {
		gwf_diag_t t;
		uint32_t v;
		int32_t d, k, i, vl;

		t = A_shift();
		v = t.vd >> 16;
		d = (int32_t)(t.vd & 0xFFFF) - GWF_DIAG_SHIFT;
		k = t.k;
		vl = s_sub_copy.seq_len[v];
		k = gwf_extend1(d, k, vl,
			s_sub_copy.graphSeq, s_sub_copy.seq_off[v],
			s_ql, s_q);
		i = k + d;

		if (k + 1 < vl && i + 1 < s_ql) {
			gwf_diag_push(B_a, &B_n, v, d-1, k+1);
			gwf_diag_push(B_a, &B_n, v, d, k+1);
			gwf_diag_push(B_a, &B_n, v, d+1, k);
		} else if (i + 1 < s_ql) {
			int32_t nv = subgfa_arc_n(&s_sub_copy, v);
			int32_t j, n_ext = 0;
			const subgfa_arc_t *av = subgfa_arc_a(&s_sub_copy, v);
			gwf_intv_t *p;
			if (s_next_intv_buf_n >= INTV_CAP) {
				fprintf(stderr, "FATAL: s_next_intv_buf overflow in tiled phase2 (n=%zu, cap=%d)\n", s_next_intv_buf_n, INTV_CAP);
				exit(1);
			}
			p = &s_next_intv_buf[s_next_intv_buf_n++];
			p->vd0 = gwf_gen_vd(v, d);
			p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) {
				uint32_t w = av[j].w;
				int32_t ol = av[j].ow;
				int absent;
				ha_put((uint32_t)w<<16 | ((i + 1) & 0xFFFF), &absent);
				if (GET_2BIT(s_q, i + 1) == GET_2BIT(s_sub_copy.graphSeq, s_sub_copy.seq_off[w] + ol)) {
					++n_ext;
					if (absent) {
						gwf_diag_t *dp = A_pushp();
						dp->vd = gwf_gen_vd(w, i + 1 - ol);
						dp->k = ol;
					}
				} else if (absent) {
					gwf_diag_push(B_a, &B_n, w, i - ol, ol);
					gwf_diag_push(B_a, &B_n, w, i + 1 - ol, ol);
				}
			}
			if (nv == 0 || n_ext != nv)
				gwf_diag_push(B_a, &B_n, v, d+1, k);
		} else if (v == GWFA_END_V && k + 1 == vl) {
			s_last_score = s;
			return 1;
		} else if (k + 1 < vl) {
			gwf_diag_push(B_a, &B_n, v, d-1, k+1);
		} else {
			int32_t nv = subgfa_arc_n(&s_sub_copy, v), j;
			const subgfa_arc_t *av = subgfa_arc_a(&s_sub_copy, v);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(B_a, &B_n, av[j].w, i - av[j].ow, av[j].ow);
		}
	}

	s_a = B_a;
	s_n_a = B_n;

	if (do_dedup)
		s_n_a = gwf_dedup(s_n_a, s_a, n_sorted);

	if (s_n_a == 0) {
		fprintf(stderr, "gwfa_extend_step_tiled: n_a==0 at s=%d\n", s);
		s_last_score = -1;
		return 1;
	}
	return 0;
}

/* ---- Split API for ISA-driven tile loop ---- */
static gwf_diag_t *s_B_a;
static int32_t s_B_n;
void gwfa_begin_step(void)
{
	s_next_intv_buf_n = 0;
	ha_clear();
	A_clear();
	s_B_a = (s_a == s_diag_a) ? s_diag_b : s_diag_a;
	s_B_n = 0;
}

/* Load one fixed-stride tile into spm.
   Returns tile_n (number of diags loaded, 0..N_TILE_DIAGS). */
static int32_t tile_load_one(int32_t cursor, int *spm)
{
	int32_t remaining = s_n_a - cursor;
	if (remaining <= 0) {
		spm[META_OFF] = 0;
		spm[META_OFF + 1] = 0;
		spm[META_OFF + 3] = 0;
		return 0;
	}
	int32_t tile_n = remaining < N_TILE_DIAGS ? remaining : N_TILE_DIAGS;
	uint64_t *dst = (uint64_t *)(spm + A_TILE_OFF);
	const uint64_t *src = (const uint64_t *)&s_a[cursor];
	for (int32_t i = 0; i < tile_n; i++)
		dst[i] = src[i]; //mvd or mvqd
	spm[META_OFF + 3] = tile_n;
	return tile_n;
}

int32_t gwfa_tile_load_one(int32_t cursor, int *spm)
{
	return tile_load_one(cursor, spm);
}

void gwfa_tile_load_seq_info(int *spm)
{
	int32_t tile_n = spm[META_OFF + 3];
	if (tile_n <= 0) { spm[META_OFF + 2] = 0; return; }
	gwf_diag_t *tile_a = (gwf_diag_t *)(spm + A_TILE_OFF);
	int32_t n_nodes = 0;
	uint32_t prev_v = UINT32_MAX;
	for (int32_t i = 0; i < tile_n; i++) {
		uint32_t v = tile_a[i].vd >> 16;
		if (v != prev_v) {
			spm[SEQ_INFO_OFF + 2*n_nodes]     = (int)s_sub_copy.seq_off[v];
			spm[SEQ_INFO_OFF + 2*n_nodes + 1] = s_sub_copy.seq_len[v];
			n_nodes++;
			prev_v = v;
		}
	}
	spm[META_OFF + 2] = n_nodes;
}

/* Like phase1_fused, but reads seq_off/seq_len from a sequential
   node_info table (pairs of int) instead of indexing by vertex ID. */
// Phase 1 for SPM-tiled execution: fused extend + next-wavefront emit.
// Extends each diagonal in `a` (Landau-Vishkin), then emits new diagonals
// into `b` for the next edit distance. Diags that reach a node or query
// boundary go into `aout` for phase 2 (cross-node propagation).
//
// For a run of consecutive diags d0..dN on the same vertex, the new
// wavefront at diagonal d uses: k_new[d] = max(k[d-1], k[d]+1, k[d+1])
// corresponding to insertion, substitution, and deletion edits.
// We compute this in a single forward pass using a sliding window:
//   ppk = k[d-2], pk = k[d-1], a[i].k = k[d]
//
// a:         input diags (sorted by vd), extended in-place
// b/bn:      output next-wavefront diags (appended)
// aout/aoutn: diags that reached a boundary (appended)
// node_info: packed [seq_offset, seq_len] per node
static void phase1_fused_spm(
	gwf_diag_t *a, int32_t n,
	gwf_diag_t *b, int32_t *bn,
	gwf_diag_t *aout, int32_t *aoutn,
	int *tile_intv, int32_t *tile_intv_n,
	const int *node_info, int32_t n_nodes,
	const uint32_t *graphSeq, const uint32_t *q, int32_t ql)
{
	int32_t i = 0, node_idx = -1;
	uint32_t prev_v = UINT32_MAX;
	while (i < n) {
		// look up node sequence offset and length
		int32_t v = a[i].vd >> 16;
		if (v != prev_v) { node_idx++; prev_v = v; }
		int32_t ts_off = node_info[2 * node_idx];
		int32_t vl = node_info[2 * node_idx + 1];
		int32_t ppk, pk, k, d;

		// --- first diagonal in this run ---
		a[i].k = gwf_extend1(
			(int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
			a[i].k, vl, graphSeq, ts_off, ql, q);
		// emit deletion: new diag at d-1 with k+1
		emit_b_tile(b, bn, tile_intv, tile_intv_n, a[i].vd - 1, a[i].k + 1, v, vl, ql);
		d = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
		if (a[i].k == vl - 1 || d + a[i].k == ql - 1)
			aout[(*aoutn)++] = a[i];
		ppk = a[i].k;
		i++;

		if (i < n && a[i].vd == a[i-1].vd + 1) {
			// --- second diagonal (have ppk, need sub between d0,d1) ---
			a[i].k = gwf_extend1(
				(int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
				a[i].k, vl, graphSeq, ts_off, ql, q);
			// emit substitution at d0: max(k[d0], k[d1]) + 1
			k = (ppk > a[i].k ? ppk : a[i].k) + 1;
			emit_b_tile(b, bn, tile_intv, tile_intv_n, a[i-1].vd, k, v, vl, ql);
			d = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
			if (a[i].k == vl - 1 || d + a[i].k == ql - 1)
				aout[(*aoutn)++] = a[i];
			pk = ppk;
			ppk = a[i].k;
			i++;

			// --- 3+ consecutive diags: full sliding window ---
			while (i < n && a[i].vd == a[i-1].vd + 1) {
				a[i].k = gwf_extend1(
					(int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT,
					a[i].k, vl, graphSeq, ts_off, ql, q);
				// emit at d-1: max(k[d-2], k[d-1]+1, k[d]+1)
				k = pk;
				if (ppk + 1 > k) k = ppk + 1;
				if (a[i].k + 1 > k) k = a[i].k + 1;
				emit_b_tile(b, bn, tile_intv, tile_intv_n, a[i-1].vd, k, v, vl, ql);
				d = (int32_t)(a[i].vd & 0xFFFF) - GWF_DIAG_SHIFT;
				if (a[i].k == vl - 1 || d + a[i].k == ql - 1)
					aout[(*aoutn)++] = a[i];
				// slide window forward
				pk = ppk;
				ppk = a[i].k;
				i++;
			}

			// emit substitution at last diagonal in the run
			k = pk > ppk + 1 ? pk : ppk + 1;
			emit_b_tile(b, bn, tile_intv, tile_intv_n, a[i-1].vd, k, v, vl, ql);
		} else {
			// single diagonal: emit substitution (only neighbor is self)
			emit_b_tile(b, bn, tile_intv, tile_intv_n, a[i-1].vd, ppk + 1, v, vl, ql);
		}

		// emit insertion: new diag at d+1 with k[last]
		emit_b_tile(b, bn, tile_intv, tile_intv_n, a[i-1].vd + 1, a[i-1].k, v, vl, ql);
	}
}

void gwfa_tile_compute(int *spm)
{
	int32_t tile_n = spm[META_OFF + 3];
	if (tile_n <= 0) return;
	int32_t n_nodes = spm[META_OFF + 2];

	int *tile_a = spm + A_TILE_OFF;
	int *tile_b = spm + B_TILE_OFF;
	int *tile_aout = spm + A_OUT_OFF;
	int *node_info = spm + SEQ_INFO_OFF;
	int *tile_intv = spm + INTV_TILE_OFF;
	const uint32_t *graphSeq = s_sub_copy.graphSeq;
	const uint32_t *q = s_q;
	int32_t ql = s_ql;
	int32_t tb_n = 0, ta_n = 0, intv_n = 0;

	int32_t gwf_extend1_local(int32_t d, int32_t k, int32_t vl,
		const uint32_t *ts, int32_t ts_off, int32_t ql_arg, const uint32_t *qs) {
		int32_t max_k = (ql_arg - d < vl ? ql_arg - d : vl) - 1;
		while (k < max_k && GET_2BIT(ts, ts_off + k + 1) == GET_2BIT(qs, d + k + 1))
			++k;
		return k;
	}

	void emit_b_tile_local(int *B_a, int32_t *B_n, int *t_intv, int32_t *t_intv_n,
		uint32_t vd, int32_t k, int32_t v, int32_t vl, int32_t ql_arg) {
		int32_t d = (int32_t)(vd & 0xFFFF) - 0x4000;
		if (d + k < ql_arg && k < vl) {
			int idx = (*B_n)++;
			B_a[2*idx] = (int)vd;
			B_a[2*idx+1] = k;
		} else if (k == vl) {
			uint32_t vd0 = ((uint32_t)v << 16) | (0x4000 + d);
			t_intv[2*(*t_intv_n)] = (int)vd0;
			t_intv[2*(*t_intv_n)+1] = (int)(vd0 + 1);
			(*t_intv_n)++;
		}
	}

	int32_t i = 0, node_idx = -1;
	uint32_t prev_v = UINT32_MAX;
	while (i < tile_n) {
		int32_t v = (uint32_t)tile_a[2*i] >> 16;
		if (v != prev_v) {
			node_idx++; prev_v = v;
		}
		int32_t ts_off = node_info[2 * node_idx];
		int32_t vl = node_info[2 * node_idx + 1];
		int32_t ppk, pk, k, d;

		tile_a[2*i+1] = gwf_extend1_local(
			(int32_t)((uint32_t)tile_a[2*i] & 0xFFFF) - GWF_DIAG_SHIFT,
			tile_a[2*i+1], vl, graphSeq, ts_off, ql, q);
		emit_b_tile_local(tile_b, &tb_n, tile_intv, &intv_n,
			(uint32_t)tile_a[2*i] - 1, tile_a[2*i+1] + 1, v, vl, ql);
		d = (int32_t)((uint32_t)tile_a[2*i] & 0xFFFF) - GWF_DIAG_SHIFT;
		if (tile_a[2*i+1] == vl - 1 || d + tile_a[2*i+1] == ql - 1) {
			tile_aout[2*ta_n] = tile_a[2*i];
			tile_aout[2*ta_n+1] = tile_a[2*i+1];
			ta_n++;
		}
		ppk = tile_a[2*i+1];
		i++;

		if (i < tile_n && (uint32_t)tile_a[2*i] == (uint32_t)tile_a[2*(i-1)] + 1) {
			tile_a[2*i+1] = gwf_extend1_local(
				(int32_t)((uint32_t)tile_a[2*i] & 0xFFFF) - GWF_DIAG_SHIFT,
				tile_a[2*i+1], vl, graphSeq, ts_off, ql, q);
			k = (ppk > tile_a[2*i+1] ? ppk : tile_a[2*i+1]) + 1;
			emit_b_tile_local(tile_b, &tb_n, tile_intv, &intv_n,
				(uint32_t)tile_a[2*(i-1)], k, v, vl, ql);
			d = (int32_t)((uint32_t)tile_a[2*i] & 0xFFFF) - GWF_DIAG_SHIFT;
			if (tile_a[2*i+1] == vl - 1 || d + tile_a[2*i+1] == ql - 1) {
				tile_aout[2*ta_n] = tile_a[2*i];
				tile_aout[2*ta_n+1] = tile_a[2*i+1];
				ta_n++;
			}
			pk = ppk;
			ppk = tile_a[2*i+1];
			i++;

			while (i < tile_n
				&& (uint32_t)tile_a[2*i] == (uint32_t)tile_a[2*(i-1)] + 1) {
				tile_a[2*i+1] = gwf_extend1_local(
					(int32_t)((uint32_t)tile_a[2*i] & 0xFFFF) - GWF_DIAG_SHIFT,
					tile_a[2*i+1], vl, graphSeq, ts_off, ql, q);
				k = pk;
				if (ppk + 1 > k) k = ppk + 1;
				if (tile_a[2*i+1] + 1 > k) k = tile_a[2*i+1] + 1;
				emit_b_tile_local(tile_b, &tb_n, tile_intv, &intv_n,
					(uint32_t)tile_a[2*(i-1)], k, v, vl, ql);
				d = (int32_t)((uint32_t)tile_a[2*i] & 0xFFFF) - GWF_DIAG_SHIFT;
				if (tile_a[2*i+1] == vl - 1 || d + tile_a[2*i+1] == ql - 1) {
					tile_aout[2*ta_n] = tile_a[2*i];
					tile_aout[2*ta_n+1] = tile_a[2*i+1];
					ta_n++;
				}
				pk = ppk;
				ppk = tile_a[2*i+1];
				i++;
			}

			k = pk > ppk + 1 ? pk : ppk + 1;
			emit_b_tile_local(tile_b, &tb_n, tile_intv, &intv_n,
				(uint32_t)tile_a[2*(i-1)], k, v, vl, ql);
		} else {
			emit_b_tile_local(tile_b, &tb_n, tile_intv, &intv_n,
				(uint32_t)tile_a[2*(i-1)], ppk + 1, v, vl, ql);
		}

		emit_b_tile_local(tile_b, &tb_n, tile_intv, &intv_n,
			(uint32_t)tile_a[2*(i-1)] + 1, tile_a[2*(i-1)+1], v, vl, ql);
	}

	spm[META_OFF] = tb_n;
	spm[META_OFF + 1] = ta_n;
	spm[META_OFF + 7] = intv_n;
}

void gwfa_B_push(uint32_t vd, int32_t k)
{
	s_B_a[s_B_n].vd = vd;
	s_B_a[s_B_n].k = k;
	s_B_n++;
}

void gwfa_tile_writeback_one(int *spm)
{
	int32_t tb_n = spm[META_OFF];
	int32_t ta_n = spm[META_OFF + 1];
	int32_t n_intv = spm[META_OFF + 7];
	gwf_diag_t *tile_b = (gwf_diag_t *)(spm + B_TILE_OFF);
	gwf_diag_t *tile_aout = (gwf_diag_t *)(spm + A_OUT_OFF);
	if (tb_n > 0) {
		memcpy(&s_B_a[s_B_n], tile_b, tb_n * sizeof(gwf_diag_t));
		s_B_n += tb_n;
	}
	for (int32_t j = 0; j < ta_n; j++)
		*A_pushp() = tile_aout[j];
	for (int32_t j = 0; j < n_intv; j++) {
		s_next_intv_buf[s_next_intv_buf_n].vd0 = (uint32_t)spm[INTV_TILE_OFF + 2*j];
		s_next_intv_buf[s_next_intv_buf_n].vd1 = (uint32_t)spm[INTV_TILE_OFF + 2*j + 1];
		s_next_intv_buf_n++;
	}
}

int gwfa_phase2(int32_t s)
{
	int do_dedup = 1;
	int32_t n_sorted = 0;

	while (A_size()) {
		gwf_diag_t t;
		uint32_t v;
		int32_t d, k, i, vl;

		t = A_shift();
		v = t.vd >> 16;
		d = (int32_t)(t.vd & 0xFFFF) - GWF_DIAG_SHIFT;
		k = t.k;
		vl = s_sub_copy.seq_len[v];
		k = gwf_extend1(d, k, vl,
			s_sub_copy.graphSeq, s_sub_copy.seq_off[v],
			s_ql, s_q);
		i = k + d;

		// Mid-vertex, mid-query: generate all 3 WFA edits (sub, ins, del)
		if (k + 1 < vl && i + 1 < s_ql) {
			gwf_diag_push(s_B_a, &s_B_n, v, d-1, k+1);
			gwf_diag_push(s_B_a, &s_B_n, v, d, k+1);
			gwf_diag_push(s_B_a, &s_B_n, v, d+1, k);
		// End of vertex, mid-query: cross arcs to successor vertices
		} else if (i + 1 < s_ql) {
			int32_t nv = subgfa_arc_n(&s_sub_copy, v);
			int32_t j, n_ext = 0;
			const subgfa_arc_t *av = subgfa_arc_a(&s_sub_copy, v);
			gwf_intv_t *p;
			if (s_next_intv_buf_n >= INTV_CAP) {
				fprintf(stderr, "FATAL: s_next_intv_buf overflow in phase2\n");
				exit(1);
			}
			// Record forbidden interval for dedup on this vertex/diagonal
			p = &s_next_intv_buf[s_next_intv_buf_n++];
			p->vd0 = gwf_gen_vd(v, d);
			p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) {
				uint32_t w = av[j].w;
				int32_t ol = av[j].ow;
				int absent;
				ha_put((uint32_t)w<<16 | ((i + 1) & 0xFFFF), &absent);
				// First char of successor matches query: extend for free
				if (GET_2BIT(s_q, i + 1) == GET_2BIT(s_sub_copy.graphSeq, s_sub_copy.seq_off[w] + ol)) {
					++n_ext;
					if (absent) {
						// Push onto active queue (same edit distance)
						gwf_diag_t *dp = A_pushp();
						dp->vd = gwf_gen_vd(w, i + 1 - ol);
						dp->k = ol;
					}
				// First char mismatches: push sub+ins edits to next wavefront
				} else if (absent) {
					gwf_diag_push(s_B_a, &s_B_n, w, i - ol, ol);
					gwf_diag_push(s_B_a, &s_B_n, w, i + 1 - ol, ol);
				}
			}
			// If not all successors matched, also push a deletion edit
			if (nv == 0 || n_ext != nv)
				gwf_diag_push(s_B_a, &s_B_n, v, d+1, k);
		// End vertex fully consumed: alignment complete, return score
		} else if (v == GWFA_END_V && k + 1 == vl) {
			s_last_score = s;
			s_a = s_B_a;
			s_n_a = s_B_n;
			return 1;
		// End of query, mid-vertex: only deletion possible (advance in ref)
		} else if (k + 1 < vl) {
			gwf_diag_push(s_B_a, &s_B_n, v, d-1, k+1);
		// End of both vertex and query, not terminal: cross to successors
		} else {
			int32_t nv = subgfa_arc_n(&s_sub_copy, v), j;
			const subgfa_arc_t *av = subgfa_arc_a(&s_sub_copy, v);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(s_B_a, &s_B_n, av[j].w, i - av[j].ow, av[j].ow);
		}
	}

	s_a = s_B_a;
	s_n_a = s_B_n;


	if (do_dedup)
		s_n_a = gwf_dedup(s_n_a, s_a, n_sorted);

	if (s_n_a == 0) {
		fprintf(stderr, "gwfa_phase2: n_a==0 at s=%d\n", s);
		s_last_score = -1;
		return 1;
	}
	return 0;
}

int32_t gwfa_get_n_a(void)
{
	return s_n_a;
}

void gwfa_debug_step(int32_t s)
{
	if (s_dbg >= 1) {
		FILE *fp = gwf_wf_debug_fp();
		fprintf(fp, "[gfa_ed_step] dist=%d, n=%d, n_intv=%zd, n_tb=0\n",
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

/* ---- MM access API ---- */

int *gwfa_get_mm(void)
{
	mm_init();
	return s_mm;
}

// Sync controller gr[] counters back to gwfa.c statics
void gwfa_sync_counters(int32_t B_n,
	uint32_t A_head, uint32_t A_tail,
	uint32_t A_count, int32_t intv_n)
{
	s_B_n = B_n;
	s_A_head = A_head;
	s_A_tail = A_tail;
	s_A_count = A_count;
	s_next_intv_buf_n = (size_t)intv_n;
}

// Return MM word offset of current s_a
int32_t gwfa_get_s_a_mm_off(void)
{
	return (int32_t)((int *)s_a - s_mm);
}

// Return MM word offset of current s_B_a
int32_t gwfa_get_s_B_a_mm_off(void)
{
	return (int32_t)((int *)s_B_a - s_mm);
}

int32_t gwfa_get_mm_A_off(void) { return MM_A_OFF; }
int32_t gwfa_get_mm_intv_off(void)
{
	return MM_NEXT_INTV_OFF;
}

int32_t gwfa_get_mm_ha_off(void) { return MM_HA_OFF; }
int32_t gwfa_get_mm_ha_dirty_off(void)
{
	return MM_HA_DIRTY_OFF;
}

// Read gwfa.c statics into caller-provided pointers
void gwfa_read_counters(int32_t *B_n,
	uint32_t *A_head, uint32_t *A_tail,
	uint32_t *A_count, int32_t *intv_n)
{
	*B_n = s_B_n;
	*A_head = s_A_head;
	*A_tail = s_A_tail;
	*A_count = s_A_count;
	*intv_n = (int32_t)s_next_intv_buf_n;
}

void gwfa_set_ha_n_dirty(uint32_t n)
{
	s_ha_n_dirty = n;
}

uint32_t gwfa_get_ha_n_dirty(void)
{
	return s_ha_n_dirty;
}

/* ---- Phase 2 tiled API ---- */

int32_t gwfa_A_size(void)
{
	return (int32_t)A_size();
}

// Pop up to P2_TILE_SIZE diags from A into spm.
// Writes (vd, k, ts_off, vl) per diag at P2_INPUT_OFF.
// Sets spm[P2_META_OFF + P2_M_TILE_N] = tile_n.
// Returns tile_n.
int32_t gwfa_phase2_tile_load(int *spm)
{
	int32_t avail = (int32_t)A_size();
	int32_t tile_n = avail < P2_TILE_SIZE ? avail : P2_TILE_SIZE;
	for (int32_t i = 0; i < tile_n; i++) {
		gwf_diag_t t = A_shift();
		uint32_t v = t.vd >> 16;
		spm[P2_INPUT_OFF + 4*i]     = (int)t.vd;
		spm[P2_INPUT_OFF + 4*i + 1] = t.k;
		spm[P2_INPUT_OFF + 4*i + 2] = (int)s_sub_copy.seq_off[v];
		spm[P2_INPUT_OFF + 4*i + 3] = s_sub_copy.seq_len[v];
	}
	spm[P2_META_OFF + P2_M_TILE_N] = tile_n;
	return tile_n;
}

// Process results from one PE's SPM after magic 13.
// Copies pushed diags to s_B_a, intvs to s_next_intv_buf.
// Processes finished_endQ0: arc traversal with ha_put.
// Processes finished_endQ1: pushes all arcs to B.
void gwfa_phase2_tile_writeback(int *spm)
{
	int32_t n_pushed = spm[P2_META_OFF + P2_M_PUSHED];
	int32_t n_intv   = spm[P2_META_OFF + P2_M_INTV];
	int32_t n_fin0   = spm[P2_META_OFF + P2_M_FIN0];
	int32_t n_fin1   = spm[P2_META_OFF + P2_M_FIN1];

	// Copy pushed diags to s_B_a
	for (int32_t j = 0; j < n_pushed; j++) {
		s_B_a[s_B_n].vd =
			(uint32_t)spm[P2_PUSHED_OFF + 2*j];
		s_B_a[s_B_n].k =
			spm[P2_PUSHED_OFF + 2*j + 1];
		s_B_n++;
	}

	// Copy intervals to s_next_intv_buf
	for (int32_t j = 0; j < n_intv; j++) {
		s_next_intv_buf[s_next_intv_buf_n].vd0 =
			(uint32_t)spm[P2_INTV_OFF + 2*j];
		s_next_intv_buf[s_next_intv_buf_n].vd1 =
			(uint32_t)spm[P2_INTV_OFF + 2*j + 1];
		s_next_intv_buf_n++;
	}

	// Process finished endQ=0: arc traversal
	for (int32_t j = 0; j < n_fin0; j++) {
		uint32_t vd = (uint32_t)spm[P2_FIN0_OFF + 2*j];
		int32_t k   = spm[P2_FIN0_OFF + 2*j + 1];
		uint32_t v  = vd >> 16;
		int32_t d   = (int32_t)(vd & 0xFFFF) - GWF_DIAG_SHIFT;
		int32_t i_val = d + k;
		int32_t nv  = subgfa_arc_n(&s_sub_copy, v);
		const subgfa_arc_t *av = subgfa_arc_a(&s_sub_copy, v);
		int32_t n_ext = 0;
		for (int32_t a = 0; a < nv; a++) {
			uint32_t w  = av[a].w;
			int32_t ol  = av[a].ow;
			int absent;
			ha_put((uint32_t)w << 16 | ((i_val + 1) & 0xFFFF), &absent);
			if (GET_2BIT(s_q, i_val + 1)
				== GET_2BIT(s_sub_copy.graphSeq, s_sub_copy.seq_off[w] + ol))
			{
				n_ext++;
				if (absent) {
					gwf_diag_t *dp = A_pushp();
					dp->vd = gwf_gen_vd(w, i_val + 1 - ol);
					dp->k = ol;
				}
			} else if (absent) {
				gwf_diag_push(s_B_a, &s_B_n, w, i_val - ol, ol);
				gwf_diag_push(s_B_a, &s_B_n, w, i_val + 1 - ol, ol);
			}
		}
		if (nv == 0 || n_ext != nv)
			gwf_diag_push(s_B_a, &s_B_n, v, d + 1, k);
	}

	// Process finished endQ=1: push all arcs to B
	for (int32_t j = 0; j < n_fin1; j++) {
		uint32_t vd = (uint32_t)spm[P2_FIN1_OFF + 2*j];
		int32_t k   = spm[P2_FIN1_OFF + 2*j + 1];
		uint32_t v  = vd >> 16;
		int32_t d   = (int32_t)(vd & 0xFFFF) - GWF_DIAG_SHIFT;
		int32_t i_val = d + k;
		int32_t nv  = subgfa_arc_n(&s_sub_copy, v);
		const subgfa_arc_t *av = subgfa_arc_a(&s_sub_copy, v);
		for (int32_t a = 0; a < nv; a++)
			gwf_diag_push(s_B_a, &s_B_n, av[a].w, i_val - av[a].ow, av[a].ow);
	}
}

// After phase 2 tile loop: dedup, set n_a.
// Returns 1 if n_a==0 (error), 0 otherwise.
int gwfa_phase2_finalize(void)
{
	s_a = s_B_a;
	s_n_a = s_B_n;
	s_n_a = gwf_dedup(s_n_a, s_a, 0);
	if (s_n_a == 0) {
		fprintf(stderr,
			"gwfa_phase2_finalize: n_a==0\n");
		s_last_score = -1;
		return 1;
	}
	return 0;
}

size_t gwfa_get_intv_n(void) { return s_intv_n; }

// Sync finalize results from controller (no dedup logic)
void gwfa_finalize_sync(int32_t n_a, size_t intv_n)
{
	s_a = s_B_a;
	s_n_a = n_a;
	s_intv_n = intv_n;
}

void gwfa_set_score(int32_t s)
{
	s_last_score = s;
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
	free(sub->arc_off);
	free(sub);
}
