# mMIG Optimization Flow - ESL Paper Artifact

Artifact for the paper:
This repository contains the pre-computed results, the figure generation script,
and the scripts to fully reproduce Table 1 from the paper.

---

## What the paper does

An mMIG is a logic graph where each internal node is either a majority gate
(MAJ) or a minority gate (MIN). This matters for AQFP superconducting circuits,
where every inverter on a fanout edge costs Josephson junctions. The paper
implements four passes layered on top of MIG optimization in
[mockturtle](https://github.com/lsils/mockturtle):

1. **Minority seeding** — seeds MAJ→MIN conversions ranked by predicted inverter
   gain across the immediate fanin/fanout neighborhood.
2. **Dual-inversion propagation** — uses MAJ/MIN complement identities to move
   inversions toward locations where a minority gate can absorb them.
3. **Cone polarity flip** — flips connected MAJ/MIN cones when the boundary
   polarity change reduces complemented fanout edges.
4. **CEC-guarded DSE** — wraps the flow in a design-space exploration loop
   that keeps only formally equivalent results (verified with ABC's `cec`).

The reported AQFP cost is the pre-mapping proxy `A = 6G + 2I`, where `G` is
gate count and `I` is complemented-edge count. This proxy counts the gate and
inverter contribution only; it does not include post-mapping splitters,
buffers, placement, or routing.

---

## Pre-computed results

`data/results.json` contains the exact numbers from Table 1. Each row is a
benchmark from the EPFL combinational suite. The comparison is *paired*: the
pure-MIG seed and the mMIG result come from the same DSE candidate, so the only
variable is the mMIG post-processing stage.

Summary across 12 benchmarks:

| Metric | Pure MIG | mMIG | Reduction |
|--------|----------|------|-----------|
| Gates | 16727 | 16285 | −2.6% |
| Inverters | 13179 | 11547 | −12.4% |
| AQFP JJ proxy | 126720 | 120804 | −4.7% |

To regenerate the figure from the pre-computed data (no build needed):

```bash
python3 scripts/gen_figure.py
# outputs fig/mmig_comparison.pdf and fig/mmig_comparison.png
```

---

## Full reproduction

### Dependencies

- CMake ≥ 3.16 and a C++17 compiler
- [ABC](https://github.com/berkeley-abc/abc) built and available as `abc` on PATH
- Python 3 with `matplotlib` and `numpy`
- EPFL combinational benchmarks (`.blif` format)

The EPFL benchmarks are freely available from the
[EPFL LSI group](https://www.epfl.ch/labs/lsi/page-102566-en-html/). Download
and extract the arithmetic and random/control suites; the 12 benchmarks used
are: adder, arbiter, bar, cavlc, ctrl, dec, i2c, int2float, max, priority,
router, voter.

### Build and run

Set `MOCKTURTLE_DIR` to the mMIG mockturtle source tree and `BENCH_DIR` to the
directory containing the EPFL BLIF files:

```bash
export MOCKTURTLE_DIR=/path/to/Mockturtle-mMIG-main
export BENCH_DIR=/path/to/epfl/blif/files
bash scripts/reproduce.sh
```

If this artifact directory is placed inside the mockturtle source tree,
`MOCKTURTLE_DIR` defaults to the parent directory and can be omitted.

The script will:
1. Build the `blif2mig_2` and `equ` binaries from the mockturtle fork.
2. Run the DSE for all 12 benchmarks (~30–60 min on a 16-core machine).
3. Extract results into `data/results.json`.
4. Regenerate the figure.

If you only want to regenerate the figure from the pre-computed table, skip the
DSE and run:

```bash
python3 scripts/gen_figure.py
```

`scripts/extract_results.py` is for full reproduction runs. It expects DSE logs
under `results/raw/<benchmark>/dse_rows.json`, which are produced by
`scripts/reproduce.sh`.

---

## Repository layout

```
data/
  results.json              — per-benchmark metrics (Table 1 numbers)
fig/
  mmig_comparison.pdf       — pre-generated comparison figure
  mmig_comparison.png       — pre-generated comparison figure
paper/
  mmig_esl.pdf              — submitted manuscript
results/
  raw/                      — per-benchmark DSE logs (populated by reproduce.sh)
scripts/
  reproduce.sh              — end-to-end reproduction script
  extract_results.py        — parse DSE logs → data/results.json
  gen_figure.py             — data/results.json → fig/mmig_comparison.*
src/
  include/mockturtle/algorithms/mmig_*.hpp  — all mMIG algorithm headers
  examples/blif2mig_2.cpp   — main binary source
  tools/mmig_dse_explore.py — design-space exploration script
  test/                     — unit tests for all mMIG passes
  configs/flows/            — flow configuration JSON files
  README.md                 — how to apply src/ to a mockturtle checkout
```

---
## Citation

The EPFL benchmark suite:

```bibtex
@inproceedings{epfl2015,
  author    = {Amarù, Luca and Gaillardon, Pierre-Emmanuel and De Micheli, Giovanni},
  title     = {The {EPFL} Combinational Benchmark Suite},
  booktitle = {Proceedings of International Workshop on Logic Synthesis (IWLS)},
  year      = {2015}
}
```
