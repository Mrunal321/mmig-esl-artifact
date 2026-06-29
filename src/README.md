# mMIG Source Files

These are the mMIG-specific source files from this work, designed to be
dropped into a [mockturtle](https://github.com/lsils/mockturtle) checkout.

## How to use

1. Clone upstream mockturtle:
   ```bash
   git clone https://github.com/lsils/mockturtle.git
   cd mockturtle
   ```

2. Copy the mMIG files over:
   ```bash
   cp -r path/to/this/src/include/mockturtle/algorithms/mmig_*.hpp \
         include/mockturtle/algorithms/
   cp path/to/this/src/include/mockturtle/algorithms/mig_resub.hpp \
         include/mockturtle/algorithms/
   cp path/to/this/src/include/mockturtle/algorithms/resubstitution.hpp \
         include/mockturtle/algorithms/
   cp path/to/this/src/examples/blif2mig_2.cpp examples/
   cp path/to/this/src/tools/mmig_dse_explore.py tools/
   ```

3. Build:
   ```bash
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release -DMOCKTURTLE_EXAMPLES=ON
   make blif2mig_2 equ -j$(nproc)
   ```

4. Run on a benchmark:
   ```bash
   ./examples/blif2mig_2 your_circuit.blif output_dir \
     --mode=area --mig-flow=dac19_compat \
     --enable-mmig --mmig-flow=compress2rs \
     --mmig-advanced --mmig-advanced-rounds=2 --mmig-cec
   ```

## What each file does

| File | Description |
|------|-------------|
| `mmig_optimizer.hpp` | Top-level orchestrator: calls seeding, rewriting, resubstitution in sequence |
| `mmig_minority_seeding.hpp` | Seeds MAJ→MIN conversions ranked by predicted inverter gain (L1/L2 lookahead) |
| `mmig_inv_propagation.hpp` | Propagates inversions through MAJ/MIN using complement identities |
| `mmig_algebraic_rewriting.hpp` | Algebraic rewriting rules specific to mixed MAJ/MIN graphs |
| `mmig_exact_rewriting.hpp` | Exact rewriting for small mMIG windows |
| `mmig_resubstitution.hpp` | Resubstitution pass adapted for mMIG inverter-aware objective |
| `mmig_cone_polarity_flip.hpp` | Flips connected MAJ/MIN cones when boundary polarity change saves inverters |
| `mmig_cec_guard.hpp` | Wraps any transform with a miter-based equivalence check |
| `mmig_cut_rewriting.hpp` | Cut-based rewriting for mMIG |
| `mmig_refactoring.hpp` | Refactoring pass for mMIG |
| `mmig_balancing.hpp` | Depth balancing for mMIG |
| `mig_resub.hpp` | Modified upstream MIG resubstitution (minor mMIG-aware changes) |
| `resubstitution.hpp` | Modified upstream resubstitution framework |
| `blif2mig_2.cpp` | Main binary: reads BLIF, runs MIG+mMIG flow, writes optimized BLIF/Verilog |
| `mmig_dse_explore.py` | Design-space exploration script: stochastic multi-restart DSE with CEC guard |
| `mmig_ablation_report.py` | Generates ablation comparison tables from DSE logs |
