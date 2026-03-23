#ifndef GWFA_H
#define GWFA_H

#include <stdint.h>

typedef struct {
	uint16_t v, w;  // packed into one 32-bit word
	int32_t ow;
} subgfa_arc_t;     // 8 bytes = 2 words

typedef struct {
	uint32_t n_vtx;
	uint32_t n_arc;
	uint32_t *graphSeq;
	uint32_t *seq_off;
	int32_t *seq_len;
	subgfa_arc_t *arc;
	uint32_t *arc_off; // CSR: n_arcs(v) = arc_off[v+1]-arc_off[v]
} subgfa_subgraph_t;

#define GWFA_START_V  0
#define GWFA_END_V    2
#define GWF_DIAG_SHIFT 0x4000

#ifdef __cplusplus
extern "C" {
#endif

int gwfa(int32_t ql, const uint32_t *q,
	subgfa_subgraph_t *sub, int32_t s_term,
	int dbg);

void gwfa_init(int32_t ql, const uint32_t *q,
	const subgfa_subgraph_t *sub, int dbg);
void gwfa_reset_step(void);
int gwfa_extend_step(int32_t s);
int gwfa_extend_step_tiled(int32_t s, int *spm);
void gwfa_begin_step(void);
int32_t gwfa_tile_load_one(int32_t cursor, int *spm);
void gwfa_tile_load_seq_info(int *spm);
void gwfa_tile_compute(int *spm);
void gwfa_B_push(uint32_t vd, int32_t k);
void gwfa_tile_writeback_one(int *spm);
int gwfa_phase2(int32_t s);
int32_t gwfa_get_n_a(void);
void gwfa_debug_step(int32_t s);
int gwfa_get_score(void);
int32_t gwfa_A_size(void);
int32_t gwfa_phase2_tile_load(int *spm);
void gwfa_phase2_tile_writeback(int *spm);
int gwfa_phase2_finalize(void);
void gwfa_set_score(int32_t s);
int *gwfa_get_mm(void);
void gwfa_sync_counters(int32_t B_n,
	uint32_t A_head, uint32_t A_tail,
	uint32_t A_count, int32_t intv_n);
int32_t gwfa_get_s_a_mm_off(void);
int32_t gwfa_get_s_B_a_mm_off(void);
int32_t gwfa_get_mm_A_off(void);
int32_t gwfa_get_mm_intv_off(void);
void gwfa_read_counters(int32_t *B_n,
	uint32_t *A_head, uint32_t *A_tail,
	uint32_t *A_count, int32_t *intv_n);

int32_t gwfa_get_mm_ha_off(void);
int32_t gwfa_get_mm_ha_dirty_off(void);
void gwfa_set_ha_n_dirty(uint32_t n);
uint32_t gwfa_get_ha_n_dirty(void);
uint32_t gwfa_ha_put(uint32_t key, int *absent);

void subgfa_subgraph_destroy(
	subgfa_subgraph_t *sub);

#ifdef __cplusplus
}
#endif

#endif
