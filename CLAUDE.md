# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GWFA (Graph Wavefront Alignment) — a C/C++ kernel that computes edit distance
via wavefront alignment on a directed acyclic graph of genomic sequences. Part
of the GenDP2 project.

## Build & Test

```bash
make                # Build (gcc-11/g++-11, -O3)
make clean          # Remove artifacts
make DBG=1          # Debug output (-DGFA_ED_DBG=4)
make GDB=1          # Debug build (-O0 -ggdb)
make ASAN=1         # AddressSanitizer build

python3 run.py 1    # Quick test (15 queries)
python3 run.py 2    # Full test (all queries)
python3 run.py 3    # Debug trace test (wfDebug output)
```

Tests rebuild from clean, run `./gwfa`, diff output against ground truth in
`../GroundTruth/`, and report kernel timing speedup.

## Architecture

**Data flow:** Text dump files → `main.cpp` loader → `gwfa()` C function →
`scores.txt`

### Key files

- **gwfa.h** — Public API: `gwfa()` entry point, `subgfa_subgraph_t` graph
  struct (vertices, arcs, sequences, offsets, index).
- **gwfa.c** — Core algorithm (631 lines). Uses static pre-allocated memory
  pools (16M diagonal buffers, 4M hash slots) to avoid malloc in hot path.
  Wavefront extends diagonals along graph sequences, generates new diagonals
  at node boundaries, dedups, and checks termination.
- **main.cpp** — C++ driver that loads datasets from `Datasets/Gwfa256/`,
  builds `subgfa_subgraph_t` per query, calls `gwfa()`, writes scores.
- **ksort.h** — Macro-based radix/heap sort (from Attractive Chaos lib).
- **kvec.h** — Simple dynamic array macros.
- **run.py** — Test runner: compiles, executes, diffs against ground truth.

### Algorithm internals (gwfa.c)

Core data structures:
- `gwf_diag_t` — diagonal state: (vertex, distance offset, k-value)
- `gwf_intv_t` — forbidden diagonal intervals for dedup

Core functions:
- `gwf_ed_extend()` — main wavefront extension loop per edit distance
- `gwf_ed_extend_batch()` — batch Landau-Vishkin extension
- `gwf_extend1()` — character-by-character match extension
- `gwf_diag_dedup()` — remove redundant diagonals via hash set
- `gwf_intv_merge_adj()` — merge overlapping forbidden intervals

Two ping-pong diagonal buffers (`s_diag_a`/`s_diag_b`) swap each iteration.
Active queue (`s_A`) uses circular buffer with shift-based indexing.

### Data format

Input lives in `Datasets/Gwfa256/` as line-per-query text files: query
sequences (`q.txt`), lengths (`ql.txt`), graph topology (`arc_v.txt`,
`arc_w.txt`, `arc_ow.txt`), vertex sequences (`graphSeq.txt`, `seq_off.txt`,
`seq_len.txt`), start/end positions, and termination steps.
