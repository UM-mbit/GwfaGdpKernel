#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include "gwfa.h"
#include "ksort.h"
#include "kvec.h"

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

#define GWF_DIAG_SHIFT 0x40000000

static inline uint64_t gwf_gen_vd(
	uint32_t v, int32_t d)
{
	return (uint64_t)v<<32
		| (GWF_DIAG_SHIFT + d);
}

/*
 * Diagonal interval
 */
typedef struct {
	uint64_t vd0, vd1;
} gwf_intv_t;

typedef kvec_t(gwf_intv_t) gwf_intv_v;

#define intvd_key(x) ((x).vd0)
KRADIX_SORT_INIT(gwf_intv, gwf_intv_t,
	intvd_key, 8)

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
	uint64_t st, en;
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
	uint64_t vd;
	int32_t k;
} gwf_diag_t;

typedef kvec_t(gwf_diag_t) gwf_diag_v;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

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
static uint64_t   s_ha_keys[HA_CAP];
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
	uint64_t key, int *absent)
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
	gwf_diag_v *a, uint32_t v,
	int32_t d, int32_t k)
{
	gwf_diag_t *p = &a->a[a->n++];
	p->vd = gwf_gen_vd(v, d), p->k = k;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(
	gwf_diag_t *p, uint32_t v,
	int32_t d, int32_t k)
{
	uint64_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->k = p->k > k ? p->k : k;
		return 0;
	}
	return 1;
}

// sort using n_sorted as split point
static void gwf_diag_sort(int32_t n_a,
	gwf_diag_t *a, int32_t n_sorted,
	gwf_diag_v *buf)
{
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	n_b = n_sorted;
	n_c = n_a - n_sorted;
	b = buf->a, c = b + n_b;
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
	gwf_diag_v *buf)
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
	gwf_diag_v sb;
	sb.a = s_sort_buf;
	sb.n = 0;
	sb.m = DIAG_CAP;
	n_a = gwf_diag_dedup(
		n_a, a, n_sorted, &sb);
	if (s_intv_n > 0)
		n_a = gwf_mixed_dedup(n_a, a,
			s_intv_n, s_intv);
	return n_a;
}

// reach the wavefront
static inline int32_t gwf_extend1(
	int32_t d, int32_t k,
	int32_t vl, const char *ts,
	int32_t ql, const char *qs)
{
	int32_t max_k =
		(ql - d < vl? ql - d : vl) - 1;
	const char *ts_ = ts + 1, *qs_ = qs + d + 1;
	uint64_t cmp = 0;
	while (k + 7 < max_k) {
		uint64_t x = *(uint64_t*)(ts_ + k);
		uint64_t y = *(uint64_t*)(qs_ + k);
		cmp = x ^ y;
		if (cmp == 0) k += 8;
		else break;
	}
	if (cmp)
		k += __builtin_ctzl(cmp) >> 3;
	else if (k + 7 >= max_k)
		while (k < max_k
			&& *(ts_ + k) == *(qs_ + k))
			++k;
	return k;
}

// Landau-Vishkin for linear sequences
static void gwf_ed_extend_batch(
	const subgfa_subgraph_t *sub,
	int32_t ql, const char *q, int32_t n,
	gwf_diag_t *a, gwf_diag_v *B)
{
	int32_t j, m;
	int32_t v = a->vd>>32;
	int32_t vl = sub->seq_len[v];
	const char *ts =
		sub->graphSeq + sub->seq_off[v];
	gwf_diag_t *b;

	// wfa_extend
	for (j = 0; j < n; ++j)
		a[j].k = gwf_extend1(
			(int32_t)a[j].vd - GWF_DIAG_SHIFT,
			a[j].k, vl, ts, ql, q);

	// wfa_next
	b = &B->a[B->n];
	b[0].vd = a[0].vd - 1;
	b[0].k = a[0].k + 1;
	b[1].vd = a[0].vd;
	b[1].k = (n == 1 || a[0].k > a[1].k
		? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k;
		k = k > a[j].k + 1
			? k : a[j].k + 1;
		k = k > a[j+1].k + 1
			? k : a[j+1].k + 1;
		b[j+1].vd = a[j].vd, b[j+1].k = k;
	}
	if (n >= 2) {
		b[n].vd = a[n-1].vd;
		b[n].k = a[n-2].k > a[n-1].k + 1
			? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].vd = a[n-1].vd + 1;
	b[n+1].k = a[n-1].k;

	// drop out-of-bound cells
	for (j = 0; j < n; ++j) {
		gwf_diag_t *p = &a[j];
		if (p->k == vl - 1
			|| (int32_t)p->vd
			- GWF_DIAG_SHIFT + p->k
			== ql - 1)
			*A_pushp() = *p;
	}
	for (j = 0, m = 0; j < n + 2; ++j) {
		gwf_diag_t *p = &b[j];
		int32_t d = (int32_t)p->vd
			- GWF_DIAG_SHIFT;
		if (d + p->k < ql && p->k < vl) {
			b[m++] = *p;
		} else if (p->k == vl) {
			if (s_tmp_n >= INTV_CAP) {
				fprintf(stderr,
					"FATAL: s_tmp overflow in "
					"gwf_ed_extend_batch "
					"(n=%zu, cap=%d)\n",
					s_tmp_n, INTV_CAP);
				exit(1);
			}
			gwf_intv_t *qi;
			qi = &s_tmp[s_tmp_n++];
			qi->vd0 = gwf_gen_vd(v, d);
			qi->vd1 = qi->vd0 + 1;
		}
	}
	B->n += m;
}

// wfa_extend and wfa_next combined
static gwf_diag_t *gwf_ed_extend(
	const subgfa_subgraph_t *sub,
	int32_t s, int32_t ql, const char *q,
	uint32_t endV,
	int32_t *n_a_, gwf_diag_t *a,
	int *terminate)
{
	int32_t i, x, n = *n_a_, do_dedup = 1;
	gwf_diag_v B;
	gwf_diag_t *b;

	s_tmp_n = 0;
	ha_clear();
	A_clear();

	/* B is the OTHER static buffer (ping-pong) */
	B.a = (a == s_diag_a) ? s_diag_b : s_diag_a;
	B.n = 0;
	B.m = DIAG_CAP;

	for (x = 0, i = 1; i <= n; ++i) {
		if (i == n || a[i].vd != a[i-1].vd + 1) {
			gwf_ed_extend_batch(sub, ql, q, i - x, &a[x], &B);
			x = i;
		}
	}
	if (A_size() == 0) do_dedup = 0;
	int32_t n_sorted = B.n;

	while (A_size()) {
		gwf_diag_t t;
		uint32_t v;
		int32_t d, k, i, vl;

		t = A_shift();
		v = t.vd >> 32;
		d = (int32_t)t.vd - GWF_DIAG_SHIFT;
		k = t.k;
		vl = sub->seq_len[v];
		k = gwf_extend1(d, k, vl,
			sub->graphSeq + sub->seq_off[v],
			ql, q);
		i = k + d;

		if (k + 1 < vl && i + 1 < ql) {
			gwf_diag_push(&B, v, d-1, k+1);
			gwf_diag_push(&B, v, d, k+1);
			gwf_diag_push(&B, v, d+1, k);
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
					(uint64_t)w<<32
					| (i + 1), &absent);
				if (q[i + 1]
					== (sub->graphSeq
					+ sub->seq_off[w])[ol])
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
					gwf_diag_push(&B,
						w, i - ol, ol);
					gwf_diag_push(&B,
						w, i + 1 - ol, ol);
				}
			}
			if (nv == 0 || n_ext != nv)
				gwf_diag_push(
					&B, v, d+1, k);
		} else if (endV == (uint32_t)-1
			|| (v == endV
			&& k + 1 == vl)) {
			*terminate = 1;
			return 0;
		} else if (k + 1 < vl) {
			gwf_diag_push(&B, v, d-1, k+1);
		} else {
			int32_t nv =
				subgfa_arc_n(sub, v), j;
			const subgfa_arc_t *av =
				subgfa_arc_a(sub, v);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(&B,
					av[j].w,
					i - av[j].ow,
					av[j].ow);
		}
	}

	*n_a_ = n = B.n, b = B.a;

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
			(int32_t)(a[i].vd0>>32),
			(int32_t)a[i].vd0
				- GWF_DIAG_SHIFT,
			(int32_t)a[i].vd1
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
			(int32_t)(a[i].vd >> 32);
		int32_t diag = (int32_t)a[i].vd
			- GWF_DIAG_SHIFT;
		fprintf(fp, "WF\t%d\t%d\t%d\n",
			nid, diag, a[i].k);
	}
	fflush(fp);
}

int gwfa(int32_t ql, const char *q,
	uint32_t startV, uint32_t endV,
	subgfa_subgraph_t *sub, int32_t s_term,
	int dbg)
{
	gwf_diag_t *a;
	int32_t n_a, s;
	int terminate = 0;

	/* Reset per-alignment counters */
	s_intv_n = 0;
	s_tmp_n = 0;
	memset(s_ha_occ, 0, sizeof(s_ha_occ));
	s_ha_n_dirty = 0;

	/* Initial wavefront */
	a = s_diag_a;
	n_a = 1;
	a[0].vd = gwf_gen_vd(startV, 0);
	a[0].k = -1;

	s = 0;
	while (n_a > 0) {
		a = gwf_ed_extend(sub, s, ql, q,
			endV, &n_a, a,
			&terminate);
		if (terminate || s >= s_term) break;
		++s;
		if (dbg >= 1) {
			FILE *fp = gwf_wf_debug_fp();
			fprintf(fp,
				"[gfa_ed_step] dist=%d, n=%d,"
				" n_intv=%zd, n_tb=0\n",
				s, n_a, s_intv_n);
			fflush(fp);
			gwf_ed_print_intv(
				s_intv_n, s_intv);
			gwf_ed_print_wf(n_a, a);
		}
	}
	return terminate ? s : -1;
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
