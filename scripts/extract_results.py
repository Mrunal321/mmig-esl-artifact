#!/usr/bin/env python3
"""Parse DSE logs and regenerate data/results.json.

Reads:  results/raw/<benchmark>/dse_rows.json
Writes: data/results.json

Run from the repo root:
    python3 scripts/extract_results.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

BENCHES = [
    "adder", "arbiter", "bar", "cavlc", "ctrl", "dec",
    "i2c", "int2float", "max", "priority", "router", "voter",
]

MIG_DSE = {
    "adder": (384, 129), "arbiter": (792, 108), "bar": (1906, 15),
    "cavlc": (374, 16),  "ctrl": (60, 8),       "dec": (304, 3),
    "i2c": (636, 16),    "int2float": (115, 9),  "max": (1939, 172),
    "priority": (337, 23), "router": (97, 13),   "voter": (3894, 32),
}

SEED_RE  = re.compile(r"seed=([^ ]+) gates=(\d+) depth=(\d+) inv=(\d+)")
MMIG_RE  = re.compile(r"mmig=([^ ]+) mmig_gates=(\d+|None) mmig_depth=(\d+|None) mmig_inv=(\d+|None)")
MIN_RE   = re.compile(r"\[mMIG\] minority nodes: \d+ -> (\d+)")
GTYPE_RE = re.compile(r"GateTypes: MAJ=\d+\s+MIN=(\d+)")


def pct(after: int, before: int) -> float:
    return 0.0 if before == 0 else (after - before) * 100.0 / before


def minority_nodes_from_log(log_path: Path) -> int:
    if not log_path.exists():
        return 0
    text = log_path.read_text(errors="ignore")
    m = MIN_RE.findall(text)
    if m:
        return int(m[-1])
    g = GTYPE_RE.findall(text)
    return int(g[-1]) if g else 0


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    raw  = repo / "results" / "raw"
    rows = []

    for bench in BENCHES:
        dse_rows = json.loads((raw / bench / "dse_rows.json").read_text())
        best_rows = [r for r in dse_rows if str(r.get("is_best_global")) == "1"]
        if not best_rows:
            raise RuntimeError(f"no is_best_global row for {bench}")
        best = best_rows[-1]
        if best.get("equiv_status") != "equivalent":
            raise RuntimeError(f"global-best row not equivalent for {bench}: {best.get('equiv_status')}")

        note  = best.get("note", "")
        seed  = SEED_RE.search(note)
        if seed:
            pure_g, pure_d, pure_i = map(int, seed.group(2, 3, 4))
        else:
            pure_g = int(best["metric_gates"])
            pure_d = int(best["metric_depth"])
            pure_i = int(best["metric_inv_edges"])

        mmig = MMIG_RE.search(note)
        if mmig and mmig.group(2) != "None":
            mmig_g, mmig_d, mmig_i = map(int, mmig.group(2, 3, 4))
        else:
            mmig_g, mmig_d, mmig_i = pure_g, pure_d, pure_i

        pure_a = 6 * pure_g + 2 * pure_i
        mmig_a = 6 * mmig_g + 2 * mmig_i
        min_nodes = minority_nodes_from_log(repo / str(best.get("log_path", "")))
        mig_g, mig_d = MIG_DSE[bench]

        rows.append([
            bench, mig_g, mig_d,
            pure_g, pure_d, pure_i, pure_a,
            mmig_g, mmig_d, mmig_i, mmig_a,
            min_nodes,
            pct(mmig_g, pure_g),
            pct(mmig_i, pure_i),
            pct(mmig_a, pure_a),
        ])

    data = {
        "description": (
            "Paired comparison: pure-MIG seed vs mMIG result from the same "
            "CEC-verified DSE row (is_best_global=1). "
            "AQFP JJ proxy: A = 6G + 2I."
        ),
        "columns": [
            "benchmark",
            "mig_dse_gates", "mig_dse_depth",
            "pure_gates", "pure_depth", "pure_inverters", "pure_aqfp_jj",
            "mmig_gates", "mmig_depth", "mmig_inverters", "mmig_aqfp_jj",
            "minority_nodes",
            "delta_gates_pct", "delta_inverters_pct", "delta_aqfp_pct",
        ],
        "rows": rows,
        "totals": {
            "pure_g":    sum(r[3]  for r in rows),
            "pure_inv":  sum(r[5]  for r in rows),
            "pure_aqfp": sum(r[6]  for r in rows),
            "mmig_g":    sum(r[7]  for r in rows),
            "mmig_inv":  sum(r[9]  for r in rows),
            "mmig_aqfp": sum(r[10] for r in rows),
        },
    }

    out = repo / "data" / "results.json"
    out.write_text(json.dumps(data, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
