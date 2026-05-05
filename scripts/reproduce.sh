#!/usr/bin/env bash
# reproduce.sh — re-run the full experiment and regenerate Table 1 + figure.
#
# Prerequisites:
#   - CMake >= 3.16, a C++17 compiler, python3 with matplotlib/numpy
#   - ABC (https://github.com/berkeley-abc/abc) built and on PATH as 'abc'
#   - EPFL benchmarks downloaded to $BENCH_DIR (see README for URL)
#
# Usage:
#   export BENCH_DIR=/path/to/epfl_benchmarks
#   bash scripts/reproduce.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MOCKTURTLE_DIR="${MOCKTURTLE_DIR:-$REPO_ROOT/..}"
BUILD_DIR="$REPO_ROOT/build"
BINARY="$BUILD_DIR/examples/blif2mig_2"
EQU_BINARY="$BUILD_DIR/examples/equ"
BENCH_DIR="${BENCH_DIR:-}"
RESULTS_ROOT="$REPO_ROOT/results"
JOBS="${JOBS:-$(nproc)}"

# ── 1. Check prerequisites ────────────────────────────────────────────────────
if [ -z "$BENCH_DIR" ]; then
  echo "ERROR: Set BENCH_DIR to the directory containing the EPFL .blif files."
  echo "  Download from: https://www.epfl.ch/labs/lsi/page-102566-en-html/"
  exit 1
fi
command -v abc >/dev/null 2>&1 || { echo "ERROR: 'abc' not found on PATH."; exit 1; }
if [ ! -f "$MOCKTURTLE_DIR/CMakeLists.txt" ] || [ ! -f "$MOCKTURTLE_DIR/tools/mmig_dse_explore.py" ]; then
  echo "ERROR: MOCKTURTLE_DIR must point to the mMIG mockturtle source tree."
  echo "  Current value: $MOCKTURTLE_DIR"
  echo "  Example: MOCKTURTLE_DIR=/path/to/Mockturtle-mMIG-main bash scripts/reproduce.sh"
  exit 1
fi

# ── 2. Build mockturtle binary ────────────────────────────────────────────────
echo "=== Building ==="
mkdir -p "$BUILD_DIR"
cmake -S "$MOCKTURTLE_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DMOCKTURTLE_EXAMPLES=ON \
  -DMOCKTURTLE_TEST=OFF \
  -DMOCKTURTLE_EXPERIMENTS=OFF \
  > "$BUILD_DIR/cmake.log" 2>&1
cmake --build "$BUILD_DIR" --target blif2mig_2 equ -j"$JOBS"
echo "  binary: $BINARY"

# ── 3. Run DSE for each benchmark ─────────────────────────────────────────────
BENCHES=(adder arbiter bar cavlc ctrl dec i2c int2float max priority router voter)

echo "=== Running DSE (this takes ~30-60 min on a 16-core machine) ==="
for bench in "${BENCHES[@]}"; do
  echo "  $bench ..."
  mkdir -p "$RESULTS_ROOT/raw"
  python3 "$MOCKTURTLE_DIR/tools/mmig_dse_explore.py" \
    --binary        "$BINARY" \
    --equ-binary    "$EQU_BINARY" \
    --abc-bin       "$(command -v abc)" \
    --bench-dir     "$BENCH_DIR" \
    --bench         "$bench" \
    --out-root      "$RESULTS_ROOT" \
    --subject       hybrid \
    --mode          area \
    --mig-flow      dac19_compat \
    --mmig-max-iters 3 \
    --mmig-advanced-rounds 2 \
    --hybrid-mmig-mig-flow identity \
    --restarts      4 \
    --steps         20 \
    --bailout-steps 8 \
    --step-timeout-s   1200 \
    --abc-timeout-s    240 \
    --equiv-timeout-s  1800 \
    --jobs             "$JOBS" \
    --inv-weight       0.5
done

# ── 4. Regenerate data/results.json ──────────────────────────────────────────
echo "=== Extracting results ==="
python3 "$REPO_ROOT/scripts/extract_results.py"

# ── 5. Regenerate figure ──────────────────────────────────────────────────────
echo "=== Generating figure ==="
python3 "$REPO_ROOT/scripts/gen_figure.py"

echo ""
echo "Done. Results in data/results.json, figure in fig/mmig_comparison.png."
