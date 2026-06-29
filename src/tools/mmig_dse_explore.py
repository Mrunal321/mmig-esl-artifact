#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import random
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from epfl_flow_report import run_case


PAPER2019_ARGS = (
    "--mmig-advanced",
    "--mmig-advanced-rounds=1",
    "--mmig-adv-no-exact",
    "--mmig-adv-no-balance",
    "--mmig-adv-no-objective-guard",
    "--mmig-adv-no-stagnation-stop",
)


MMIG_DEFAULT_COMPRESSORS = "paper2019,compress2rs_dc_reseed,epfl_dc_reseed,rr_rank_norm"
MIG_DEFAULT_COMPRESSORS = (
    "mig_paper2019,mig_compress2rs,mig_dac19_default,mig_dac19_area,"
    "mig_dac19_area_zg,mig_dac19_area_exact,mig_dac19_area_dc,"
    "mig_rs10,mig_rs12,mig_rs12_zg,mig_legacy"
)
HYBRID_DEFAULT_COMPRESSORS = (
    "mig_paper2019+skip_mmig,mig_rs10+skip_mmig,mig_rs12_zg+skip_mmig,"
    "mig_dac19_area_exact+skip_mmig,mig_legacy+skip_mmig,"
    "mig_paper2019+paper2019,mig_paper2019+compress2rs_dc_reseed,"
    "mig_rs10+compress2rs_dc_reseed,mig_rs12_zg+paper2019,"
    "mig_rs12_zg+compress2rs_dc_reseed,mig_dac19_area_exact+paper2019,"
    "mig_dac19_area_exact+epfl_dc_reseed,mig_dac19_area_dc+epfl_dc_reseed,"
    "mig_legacy+paper2019,mig_legacy+compress2rs_dc_reseed"
)


def default_compressors(subject: str) -> str:
    if subject == "mig":
        return MIG_DEFAULT_COMPRESSORS
    if subject == "hybrid":
        return HYBRID_DEFAULT_COMPRESSORS
    return MMIG_DEFAULT_COMPRESSORS


def default_baseline_compressor(subject: str) -> str:
    if subject == "mig":
        return "mig_paper2019"
    if subject == "hybrid":
        return "mig_paper2019+skip_mmig"
    return "paper2019"


def mmig_compressor_catalog(rounds: int) -> Dict[str, Tuple[str, Tuple[str, ...]]]:
    adv = ("--mmig-advanced", f"--mmig-advanced-rounds={rounds}")
    return {
        "native": ("round_robin", ()),
        "rr_rank_norm": ("round_robin", adv),
        "paper2019": ("paper2019", PAPER2019_ARGS),
        "paper2019_guarded": ("paper2019", ("--mmig-advanced", "--mmig-advanced-rounds=1", "--mmig-adv-no-exact", "--mmig-adv-no-balance")),
        "compress2rs": ("compress2rs", adv),
        "compress2rs_dc_reseed": ("compress2rs", adv + ("--mmig-dont-cares", "--mmig-interleaved-seeding")),
        "epfl_dc_reseed": ("epfl", adv + ("--mmig-dont-cares", "--mmig-interleaved-seeding")),
    }


def mig_compressor_catalog() -> Dict[str, Tuple[str, Tuple[str, ...]]]:
    # These are pure-MIG compression variants.  They intentionally do not pass
    # --enable-mmig; they only select the majority flow and majority-flow flags.
    return {
        "mig_standard": ("standard", ()),
        "mig_legacy": ("legacy", ()),
        "mig_legacy_zg": ("legacy", ("--mig-allow-zero-gain",)),
        "mig_dac19_default": ("dac19_default", ()),
        "mig_dac19_default_zg": ("dac19_default", ("--mig-allow-zero-gain",)),
        "mig_dac19_area": ("dac19_area", ()),
        "mig_dac19_area_zg": ("dac19_area", ("--mig-allow-zero-gain",)),
        "mig_dac19_area_exact": ("dac19_area", ("--mig-allow-zero-gain", "--mig-enable-exact")),
        "mig_dac19_area_bal": ("dac19_area", ("--mig-allow-zero-gain", "--mig-enable-balance")),
        "mig_dac19_area_dc": ("dac19_area", ("--mig-allow-zero-gain", "--resub-dc", "--resub-dc-window=8")),
        "mig_compress2rs": ("compress2rs", ()),
        "mig_compress2rs_exact": ("compress2rs", ("--mig-enable-exact",)),
        "mig_compress2rs_bal": ("compress2rs", ("--mig-enable-balance",)),
        "mig_paper2019": ("dac19_compat", ()),
        "mig_rs10": ("dac19_default", ("--rs-cut=10", "--rs-insert=2", "--rw-cut-limit=16")),
        "mig_rs12": ("dac19_default", ("--rs-cut=12", "--rs-insert=2", "--rw-cut-limit=16")),
        "mig_rs12_zg": (
            "dac19_default",
            ("--rs-cut=12", "--rs-insert=2", "--rw-cut-limit=16", "--mig-allow-zero-gain"),
        ),
        "mig_rs16_zg": (
            "dac19_default",
            ("--rs-cut=16", "--rs-insert=2", "--rw-cut-limit=16", "--mig-allow-zero-gain"),
        ),
    }


def compressor_catalog(subject: str, rounds: int) -> Dict[str, Tuple[str, Tuple[str, ...]]]:
    if subject == "mig":
        return mig_compressor_catalog()
    if subject == "hybrid":
        hybrid_mmig = dict(mmig_compressor_catalog(rounds))
        hybrid_mmig["skip_mmig"] = ("skip_mmig", ())
        return {
            f"{mig_name}+{mmig_name}": ("hybrid", ())
            for mig_name in mig_compressor_catalog()
            for mmig_name in hybrid_mmig
        }
    return mmig_compressor_catalog(rounds)


ABC_DECOMPRESS_SCRIPTS = {
    "abc_resyn2": (
        "read_blif {in}; strash; balance; rewrite; refactor; balance; rewrite; rewrite -z; "
        "balance; refactor -z; rewrite -z; balance; write_blif {out}"
    ),
    "abc_resyn2rs": (
        "read_blif {in}; strash; balance; resub; resub -K 6; balance; resub -z; "
        "resub -z -K 6; balance; resub -z -K 5; balance; write_blif {out}"
    ),
    "abc_if3": "read_blif {in}; strash; if -K 3; write_blif {out}",
    "abc_if4": "read_blif {in}; strash; if -K 4; write_blif {out}",
    "abc_if5": "read_blif {in}; strash; if -K 5; write_blif {out}",
    "abc_if6": "read_blif {in}; strash; if -K 6; write_blif {out}",
    "abc_if6_mfs": "read_blif {in}; strash; if -K 6; mfs2 -W 4 -M 500; write_blif {out}",
    "abc_dch_if3": "read_blif {in}; strash; dch; if -K 3; write_blif {out}",
    "abc_dch_if4": "read_blif {in}; strash; dch; if -K 4; write_blif {out}",
    "abc_dch_if5": "read_blif {in}; strash; dch; if -K 5; write_blif {out}",
    "abc_dch_if6": "read_blif {in}; strash; dch; if -K 6; write_blif {out}",
    "abc_dch_if6_mfs": "read_blif {in}; strash; dch; if -K 6; mfs2 -W 4 -M 500; write_blif {out}",
    "abc_dc2": "read_blif {in}; strash; dc2; write_blif {out}",
    "abc_fraig": "read_blif {in}; strash; fraig; write_blif {out}",
}


FIELDNAMES = [
    "benchmark",
    "subject",
    "restart",
    "step",
    "decompressor",
    "compressor",
    "status",
    "return_code",
    "equiv_status",
    "metric_gates",
    "metric_depth",
    "metric_inv_edges",
    "wall_s",
    "compress_wall_s",
    "decompress_wall_s",
    "is_best_global",
    "is_best_restart",
    "out_base",
    "candidate_blif",
    "intermediate_blif",
    "log_path",
    "command",
    "equiv_output",
    "note",
]


def stable_seed(base_seed: int, name: str) -> int:
    digest = hashlib.sha256(name.encode("utf-8")).hexdigest()
    return base_seed + int(digest[:8], 16)


def parse_list(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def resolve_benches(bench_dir: Path, selectors: Sequence[str]) -> List[Path]:
    by_name: Dict[str, Path] = {}
    for bench in sorted(bench_dir.rglob("*.blif")):
        current = by_name.get(bench.name)
        if current is None:
            by_name[bench.name] = bench
            continue
        current_rank = (len(current.relative_to(bench_dir).parts), str(current))
        next_rank = (len(bench.relative_to(bench_dir).parts), str(bench))
        if next_rank < current_rank:
            by_name[bench.name] = bench
    benches = sorted(by_name.values(), key=lambda p: p.name)
    if not selectors:
        return benches

    available = {b.name: b for b in benches}
    selected: List[Path] = []
    seen = set()
    for item in selectors:
        base = Path(item).name
        candidates = [base] if base.endswith(".blif") else [f"{base}.blif", f"{base}_rebuilt.blif"]
        match = next((name for name in candidates if name in available), None)
        if match is None:
            print(f"[warn] benchmark selector did not match: {item}", file=sys.stderr)
            continue
        if match not in seen:
            selected.append(available[match])
            seen.add(match)
    return selected


# Inverter weight applied to the primary cost. Default 0.0 preserves the
# original strict lex order (gates, depth, inv) where inverters only break
# ties. Set --inv-weight > 0 to fold inverters into the primary key so the
# DSE can surface mMIG-favored solutions that trade a few extra gates for
# substantially fewer inverters. The mMIG-unique wins live exactly here.
_INV_WEIGHT: float = 0.0


def set_inv_weight(weight: float) -> None:
    global _INV_WEIGHT
    _INV_WEIGHT = max(0.0, float(weight))


def _weighted_primary(gates: int, inv: int) -> float:
    if _INV_WEIGHT <= 0.0:
        return float(gates)
    return float(gates) + _INV_WEIGHT * float(inv)


def objective(row: Dict[str, object], mode: str) -> Tuple[float, int, int]:
    gates = int(row.get("metric_gates") or 10**12)
    depth = int(row.get("metric_depth") or 10**12)
    inv = int(row.get("metric_inv_edges") or 10**12)
    if mode == "depth":
        return (float(depth), _weighted_primary(gates, inv), inv)
    return (_weighted_primary(gates, inv), float(depth), inv)


def metric_objective(metrics: Dict[str, object], mode: str) -> Tuple[float, int, int]:
    gates = int(metrics.get("gates") or 10**12)
    depth = int(metrics.get("depth") or 10**12)
    inv = int(metrics.get("inv_edges_total") or 10**12)
    if mode == "depth":
        return (float(depth), _weighted_primary(gates, inv), inv)
    return (_weighted_primary(gates, inv), float(depth), inv)


def better(lhs: Optional[Dict[str, object]], rhs: Optional[Dict[str, object]], mode: str) -> bool:
    if lhs is None:
        return False
    if rhs is None:
        return True
    return objective(lhs, mode) < objective(rhs, mode)


def geometric_mean(values: Iterable[float]) -> Optional[float]:
    vals = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not vals:
        return None
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


def run_equivalence(equ_binary: Path, original: Path, revised: Path, timeout_s: int) -> Tuple[str, str]:
    if not equ_binary.exists():
        return "skipped", f"equivalence binary not found: {equ_binary}"
    if not revised.exists():
        return "missing_output", f"candidate BLIF not found: {revised}"
    try:
        completed = subprocess.run(
            [str(equ_binary), str(original), str(revised)],
            capture_output=True,
            text=True,
            timeout=timeout_s if timeout_s > 0 else None,
        )
    except subprocess.TimeoutExpired:
        return "timeout", f"equivalence exceeded {timeout_s} seconds"

    output = ((completed.stdout or "") + ("\n" + completed.stderr if completed.stderr else "")).strip()
    if completed.returncode == 0 and "EQUIVALENT" in output:
        return "equivalent", output
    if completed.returncode == 2:
        return "not_equivalent", output
    return f"error_{completed.returncode}", output


def run_abc_decompress(
    abc_bin: str,
    script_name: str,
    input_blif: Path,
    output_blif: Path,
    timeout_s: int,
) -> Tuple[bool, float, str]:
    script_template = ABC_DECOMPRESS_SCRIPTS[script_name]
    script = script_template.format_map({"in": str(input_blif), "out": str(output_blif)})
    t0 = time.time()
    try:
        completed = subprocess.run(
            [abc_bin, "-c", script],
            capture_output=True,
            text=True,
            timeout=timeout_s if timeout_s > 0 else None,
        )
    except subprocess.TimeoutExpired:
        return False, time.time() - t0, f"abc timeout after {timeout_s}s"
    wall = time.time() - t0
    output = ((completed.stdout or "") + ("\n" + completed.stderr if completed.stderr else "")).strip()
    if completed.returncode != 0:
        return False, wall, f"abc rc={completed.returncode}: {output[-400:]}"
    if not output_blif.exists():
        return False, wall, "abc did not produce output BLIF"
    return True, wall, output[-400:]


def decompress(
    name: str,
    current_blif: Path,
    out_blif: Path,
    abc_bin: str,
    abc_timeout_s: int,
) -> Tuple[Path, float, str]:
    if name == "identity":
        return current_blif, 0.0, "identity"
    if name not in ABC_DECOMPRESS_SCRIPTS:
        return current_blif, 0.0, f"unknown decompressor {name}; used identity"
    if not shutil.which(abc_bin) and not Path(abc_bin).exists():
        return current_blif, 0.0, f"abc binary not found: {abc_bin}; used identity"
    ok, wall, note = run_abc_decompress(abc_bin, name, current_blif, out_blif, abc_timeout_s)
    if not ok:
        return current_blif, wall, f"{name} failed; used identity; {note}"
    return out_blif, wall, note


def run_compressor(
    binary: Path,
    input_blif: Path,
    out_base: Path,
    mode: str,
    mig_flow: str,
    subject: str,
    compressor: str,
    compressors: Dict[str, Tuple[str, Tuple[str, ...]]],
    mmig_max_iters: int,
    mmig_stage: str,
    mmig_advanced_rounds: int,
    hybrid_mmig_mig_flow: str,
    hybrid_import_lut: bool,
    hybrid_stage_guard: bool,
    timeout_s: int,
) -> Dict[str, object]:
    if subject == "hybrid":
        mig_name, mmig_name = compressor.split("+", 1)
        mig_catalog = mig_compressor_catalog()
        mmig_catalog = mmig_compressor_catalog(mmig_advanced_rounds)
        mig_stage_base = Path(str(out_base) + "_migseed")
        mig_flow_name, mig_extra_args = mig_catalog[mig_name]
        mig_result = run_case(
            binary=binary,
            bench=input_blif,
            out_base=mig_stage_base,
            mode=mode,
            mig_flow=mig_flow_name,
            use_mmig=False,
            mmig_flow="round_robin",
            mmig_max_iters=mmig_max_iters,
            mmig_stage=mmig_stage,
            mmig_cec=True,
            use_zero_gain=False,
            extra_args=list(mig_extra_args),
            timeout_s=timeout_s,
        )
        mig_seed_blif = Path(str(mig_stage_base) + "_maj_opt.blif")
        if int(mig_result["return_code"]) != 0 or not mig_seed_blif.exists():
            stdout = "\n".join(
                [
                    "=== HYBRID MIG SEED STAGE ===",
                    str(mig_result["stdout"]),
                    f"[hybrid] MIG seed output missing or failed: {mig_seed_blif}",
                ]
            )
            return {
                **mig_result,
                "stdout": stdout,
                "command": str(mig_result["command"]),
                "mmig_metrics": {},
                "hybrid_note": f"seed stage failed or missing output {mig_seed_blif}",
            }

        hybrid_best_blif = Path(str(out_base) + "_hybrid_best.blif")
        hybrid_best_verilog = Path(str(out_base) + "_hybrid_best.v")
        mig_seed_verilog = Path(str(mig_stage_base) + "_maj_opt.v")
        seed = mig_result["mig_metrics"]
        if mmig_name == "skip_mmig":
            shutil.copy2(mig_seed_blif, hybrid_best_blif)
            if mig_seed_verilog.exists():
                shutil.copy2(mig_seed_verilog, hybrid_best_verilog)
            seed_note = (
                f"[hybrid] selected=mig_seed seed={mig_name} gates={seed.get('gates')} depth={seed.get('depth')} "
                f"inv={seed.get('inv_edges_total')} blif={mig_seed_blif}"
            )
            return {
                **mig_result,
                "stdout": "\n".join(
                    [
                        "=== HYBRID MIG SEED STAGE ===",
                        str(mig_result["stdout"]),
                        seed_note,
                    ]
                ),
                "command": str(mig_result["command"]),
                "mmig_metrics": seed,
                "hybrid_seed_blif": str(mig_seed_blif),
                "hybrid_seed_metrics": seed,
                "hybrid_selected_stage": "mig_seed",
                "hybrid_note": seed_note,
            }

        mmig_flow_name, mmig_extra_args = mmig_catalog[mmig_name]
        final_extra_args = list(mmig_extra_args)
        if hybrid_import_lut:
            final_extra_args.append("--import-lut")
        mmig_result = run_case(
            binary=binary,
            bench=mig_seed_blif,
            out_base=out_base,
            mode=mode,
            mig_flow=hybrid_mmig_mig_flow,
            use_mmig=True,
            mmig_flow=mmig_flow_name,
            mmig_max_iters=mmig_max_iters,
            mmig_stage=mmig_stage,
            mmig_cec=True,
            use_zero_gain=False,
            extra_args=final_extra_args,
            timeout_s=timeout_s,
        )
        mmig_blif = Path(str(out_base) + "_mmig_opt.blif")
        mmig_verilog = Path(str(out_base) + "_mmig_opt.v")
        seed_obj = metric_objective(seed, mode)
        mmig_obj = metric_objective(mmig_result["mmig_metrics"], mode)
        keep_seed = hybrid_stage_guard and seed_obj < mmig_obj
        if keep_seed:
            shutil.copy2(mig_seed_blif, hybrid_best_blif)
            if mig_seed_verilog.exists():
                shutil.copy2(mig_seed_verilog, hybrid_best_verilog)
            selected_stage = "mig_seed"
            selected_metrics = seed
        else:
            if mmig_blif.exists():
                shutil.copy2(mmig_blif, hybrid_best_blif)
            if mmig_verilog.exists():
                shutil.copy2(mmig_verilog, hybrid_best_verilog)
            selected_stage = "mmig"
            selected_metrics = mmig_result["mmig_metrics"]
        seed_note = (
            f"[hybrid] selected={selected_stage} seed={mig_name} gates={seed.get('gates')} depth={seed.get('depth')} "
            f"inv={seed.get('inv_edges_total')} mmig={mmig_name} "
            f"mmig_gates={mmig_result['mmig_metrics'].get('gates')} mmig_depth={mmig_result['mmig_metrics'].get('depth')} "
            f"mmig_inv={mmig_result['mmig_metrics'].get('inv_edges_total')} seed_blif={mig_seed_blif}"
        )
        return {
            **mmig_result,
            "return_code": 0 if keep_seed else mmig_result["return_code"],
            "stdout": "\n".join(
                [
                    "=== HYBRID MIG SEED STAGE ===",
                    str(mig_result["stdout"]),
                    "=== HYBRID MMIG STAGE ===",
                    str(mmig_result["stdout"]),
                    seed_note,
                ]
            ),
            "wall_s": float(mig_result["wall_s"]) + float(mmig_result["wall_s"]),
            "command": f"{mig_result['command']} && {mmig_result['command']}",
            "mig_time_s": None,
            "mmig_time_s": mmig_result["mmig_time_s"],
            "mmig_metrics": selected_metrics,
            "hybrid_seed_blif": str(mig_seed_blif),
            "hybrid_seed_metrics": seed,
            "hybrid_selected_stage": selected_stage,
            "hybrid_note": seed_note,
        }

    flow, extra_args = compressors[compressor]
    use_mmig = subject == "mmig"
    return run_case(
        binary=binary,
        bench=input_blif,
        out_base=out_base,
        mode=mode,
        mig_flow=mig_flow if use_mmig else flow,
        use_mmig=use_mmig,
        mmig_flow=flow if use_mmig else "round_robin",
        mmig_max_iters=mmig_max_iters,
        mmig_stage=mmig_stage,
        mmig_cec=True,
        use_zero_gain=False,
        extra_args=list(extra_args),
        timeout_s=timeout_s,
    )


def make_row(
    benchmark: str,
    subject: str,
    restart: int,
    step: int,
    decompressor: str,
    compressor: str,
    status: str,
    result: Optional[Dict[str, object]],
    equiv_status: str,
    equiv_output: str,
    out_base: Path,
    candidate_blif: Path,
    intermediate_blif: Path,
    log_path: Path,
    decompress_wall_s: float,
    note: str,
) -> Dict[str, object]:
    metrics_key = "mmig_metrics" if subject in ("mmig", "hybrid") else "mig_metrics"
    metrics = result[metrics_key] if result is not None else {}
    row_note = note
    if result is not None and result.get("hybrid_note"):
        row_note = f"{note}; {result['hybrid_note']}"
    return {
        "benchmark": benchmark,
        "subject": subject,
        "restart": restart,
        "step": step,
        "decompressor": decompressor,
        "compressor": compressor,
        "status": status,
        "return_code": result["return_code"] if result is not None else "",
        "equiv_status": equiv_status,
        "metric_gates": metrics.get("gates", ""),
        "metric_depth": metrics.get("depth", ""),
        "metric_inv_edges": metrics.get("inv_edges_total", ""),
        "wall_s": f"{(float(result['wall_s']) if result is not None else 0.0) + decompress_wall_s:.3f}",
        "compress_wall_s": f"{float(result['wall_s']):.3f}" if result is not None else "",
        "decompress_wall_s": f"{decompress_wall_s:.3f}",
        "is_best_global": "0",
        "is_best_restart": "0",
        "out_base": str(out_base),
        "candidate_blif": str(candidate_blif),
        "intermediate_blif": str(intermediate_blif),
        "log_path": str(log_path),
        "command": result["command"] if result is not None else "",
        "equiv_output": equiv_output,
        "note": row_note,
    }


def valid_candidate(row: Dict[str, object]) -> bool:
    return (
        str(row.get("status")) == "ok"
        and str(row.get("return_code")) == "0"
        and row.get("equiv_status") in ("equivalent", "skipped")
        and row.get("metric_gates") not in ("", None)
        and row.get("metric_depth") not in ("", None)
    )


def candidate_suffix(subject: str) -> str:
    if subject == "hybrid":
        return "_hybrid_best.blif"
    return "_mmig_opt.blif" if subject == "mmig" else "_maj_opt.blif"


def run_baseline(
    bench: Path,
    work_dir: Path,
    args: argparse.Namespace,
    compressors: Dict[str, Tuple[str, Tuple[str, ...]]],
) -> Optional[Dict[str, object]]:
    if args.no_baseline:
        return None
    compressor = args.baseline_compressor
    out_base = work_dir / f"{bench.stem}_baseline_{compressor}"
    log_path = work_dir / f"{bench.stem}_baseline_{compressor}.log"
    result = run_compressor(
        Path(args.binary),
        bench,
        out_base,
        args.mode,
        args.mig_flow,
        args.subject,
        compressor,
        compressors,
        max(1, args.mmig_max_iters),
        args.mmig_stage,
        max(1, args.mmig_advanced_rounds),
        args.hybrid_mmig_mig_flow,
        not args.hybrid_no_import_lut,
        not args.hybrid_no_stage_guard,
        max(0, args.step_timeout_s),
    )
    log_path.write_text(result["stdout"], encoding="utf-8", errors="ignore")
    candidate_blif = Path(str(out_base) + candidate_suffix(args.subject))
    equiv_status, equiv_output = ("skipped", "")
    if not args.no_equivalence and int(result["return_code"]) == 0:
        equiv_status, equiv_output = run_equivalence(Path(args.equ_binary), bench, candidate_blif, max(0, args.equiv_timeout_s))
    row = make_row(
        bench.name,
        args.subject,
        -1,
        0,
        "baseline",
        compressor,
        "ok" if int(result["return_code"]) == 0 and equiv_status in ("equivalent", "skipped") else "failed",
        result,
        equiv_status,
        equiv_output,
        out_base,
        candidate_blif,
        bench,
        log_path,
        0.0,
        "baseline",
    )
    return row if valid_candidate(row) else None


def run_benchmark(task: Tuple[Path, argparse.Namespace]) -> Dict[str, object]:
    bench, args = task
    compressors = compressor_catalog(args.subject, max(1, args.mmig_advanced_rounds))
    selected_compressors = parse_list(args.compressors)
    selected_decompressors = parse_list(args.decompressors)
    for name in selected_compressors:
        if name not in compressors:
            raise ValueError(f"unknown compressor '{name}'")
    for name in selected_decompressors:
        if name != "identity" and name not in ABC_DECOMPRESS_SCRIPTS:
            raise ValueError(f"unknown decompressor '{name}'")

    bench_dir = Path(args.out_root) / "raw" / bench.stem
    bench_dir.mkdir(parents=True, exist_ok=True)
    rows_json = bench_dir / "dse_rows.json"
    best_json = bench_dir / "dse_best.json"
    if args.resume and rows_json.exists() and best_json.exists():
        return {
            "benchmark": bench.name,
            "rows": json.loads(rows_json.read_text(encoding="utf-8")),
            "best": json.loads(best_json.read_text(encoding="utf-8")),
            "resumed": True,
        }

    rng = random.Random(stable_seed(args.seed, bench.name))
    rows: List[Dict[str, object]] = []
    best_global: Optional[Dict[str, object]] = run_baseline(bench, bench_dir, args, compressors)
    if best_global is not None:
        best_global["is_best_global"] = "1"
        best_global["is_best_restart"] = "1"
        rows.append(best_global)

    for restart in range(max(1, args.restarts)):
        restart_start = time.time()
        current_blif = bench
        best_restart: Optional[Dict[str, object]] = None
        stale_steps = 0
        for step in range(1, max(1, args.steps) + 1):
            if args.restart_timeout_s > 0 and time.time() - restart_start > args.restart_timeout_s:
                break
            if stale_steps >= args.bailout_steps:
                break

            decompressor = rng.choice(selected_decompressors)
            compressor = rng.choice(selected_compressors)
            step_tag = f"r{restart:02d}_s{step:03d}_{decompressor}_{compressor}"
            intermediate_blif = bench_dir / f"{bench.stem}_{step_tag}_decomp.blif"
            out_base = bench_dir / f"{bench.stem}_{step_tag}"
            log_path = bench_dir / f"{bench.stem}_{step_tag}.log"

            decomp_blif, decomp_wall, decomp_note = decompress(
                decompressor,
                current_blif,
                intermediate_blif,
                args.abc_bin,
                max(0, args.abc_timeout_s),
            )
            result = run_compressor(
                Path(args.binary),
                decomp_blif,
                out_base,
                args.mode,
                args.mig_flow,
                args.subject,
                compressor,
                compressors,
                max(1, args.mmig_max_iters),
                args.mmig_stage,
                max(1, args.mmig_advanced_rounds),
                args.hybrid_mmig_mig_flow,
                not args.hybrid_no_import_lut,
                not args.hybrid_no_stage_guard,
                max(0, args.step_timeout_s),
            )
            log_path.write_text(result["stdout"], encoding="utf-8", errors="ignore")
            candidate_blif = Path(str(out_base) + candidate_suffix(args.subject))

            equiv_status, equiv_output = ("skipped", "")
            if not args.no_equivalence and int(result["return_code"]) == 0:
                equiv_status, equiv_output = run_equivalence(Path(args.equ_binary), bench, candidate_blif, max(0, args.equiv_timeout_s))

            status = "ok" if int(result["return_code"]) == 0 and equiv_status in ("equivalent", "skipped") else "failed"
            row = make_row(
                bench.name,
                args.subject,
                restart,
                step,
                decompressor,
                compressor,
                status,
                result,
                equiv_status,
                equiv_output,
                out_base,
                candidate_blif,
                decomp_blif,
                log_path,
                decomp_wall,
                decomp_note,
            )

            improved_restart = valid_candidate(row) and better(row, best_restart, args.mode)
            improved_global = valid_candidate(row) and better(row, best_global, args.mode)
            if improved_restart:
                if best_restart is not None:
                    best_restart["is_best_restart"] = "0"
                best_restart = row
                row["is_best_restart"] = "1"
                stale_steps = 0
            else:
                stale_steps += 1
            if improved_global:
                if best_global is not None:
                    best_global["is_best_global"] = "0"
                best_global = row
                row["is_best_global"] = "1"

            rows.append(row)
            if valid_candidate(row):
                current_blif = candidate_blif

    if best_global is None:
        best_global = {
            "benchmark": bench.name,
            "subject": args.subject,
            "status": "failed",
            "metric_gates": "",
            "metric_depth": "",
            "metric_inv_edges": "",
            "equiv_status": "",
            "out_base": "",
            "candidate_blif": "",
            "log_path": "",
        }

    rows_json.write_text(json.dumps(rows, indent=2, sort_keys=True), encoding="utf-8")
    best_json.write_text(json.dumps(best_global, indent=2, sort_keys=True), encoding="utf-8")
    return {"benchmark": bench.name, "rows": rows, "best": best_global, "resumed": False}


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def copy_best_artifacts(selected_dir: Path, best_rows: List[Dict[str, object]]) -> None:
    selected_dir.mkdir(parents=True, exist_ok=True)
    for row in best_rows:
        if not valid_candidate(row):
            continue
        out_base = Path(str(row["out_base"]))
        dst_base = selected_dir / f"{Path(str(row['benchmark'])).stem}_{row['compressor']}_{row['decompressor']}"
        for suffix in ("_hybrid_best.blif", "_hybrid_best.v", "_mmig_opt.blif", "_mmig_opt.v", "_maj_opt.blif", "_maj_opt.v"):
            src = Path(str(out_base) + suffix)
            if src.exists():
                shutil.copy2(src, Path(str(dst_base) + suffix))
        log_path = Path(str(row["log_path"]))
        if log_path.exists():
            shutil.copy2(log_path, Path(str(dst_base) + ".log"))


def write_summary(path: Path, rows: List[Dict[str, object]], best_rows: List[Dict[str, object]], args: argparse.Namespace, csv_path: Path, best_csv: Path, selected_dir: Path) -> None:
    valid_best = [row for row in best_rows if valid_candidate(row)]
    baseline_rows = [row for row in rows if int(row.get("restart", 0)) == -1 and valid_candidate(row)]
    baseline_by_bench = {str(row["benchmark"]): row for row in baseline_rows}
    pairs = [(row, baseline_by_bench[str(row["benchmark"])]) for row in valid_best if str(row["benchmark"]) in baseline_by_bench]
    improved = same = worse = 0
    for row, base in pairs:
        if objective(row, args.mode) < objective(base, args.mode):
            improved += 1
        elif objective(row, args.mode) == objective(base, args.mode):
            same += 1
        else:
            worse += 1

    def ratio(row: Dict[str, object], base: Dict[str, object], key: str) -> Optional[float]:
        try:
            rv = float(row[key])
            bv = float(base[key])
        except (KeyError, TypeError, ValueError):
            return None
        return rv / bv if bv > 0 else None

    gate_geo = geometric_mean(r for r in (ratio(row, base, "metric_gates") for row, base in pairs) if r is not None)
    depth_geo = geometric_mean(r for r in (ratio(row, base, "metric_depth") for row, base in pairs) if r is not None)
    inv_geo = geometric_mean(r for r in (ratio(row, base, "metric_inv_edges") for row, base in pairs) if r is not None)
    wall_geo = geometric_mean(r for r in (ratio(row, base, "wall_s") for row, base in pairs) if r is not None)

    def fmt(value: Optional[float]) -> str:
        return "n/a" if value is None else f"{value:.4f}"

    mig_flow_line = (
        f"- MIG pre-flow: `{args.mig_flow}`"
        if args.subject == "mmig"
        else (
            f"- MIG pre-flow: staged hybrid; seed compressor selects MIG flow; mMIG stage MIG flow `{args.hybrid_mmig_mig_flow}`"
            if args.subject == "hybrid"
            else "- MIG pre-flow: `n/a (compressor selects MIG flow)`"
        )
    )
    lines = [
        f"# {args.subject.upper()} DSE Summary",
        "",
        f"- Subject: `{args.subject}`",
        f"- Mode: `{args.mode}`",
        mig_flow_line,
    ]
    if args.subject == "hybrid":
        lines.append(f"- Hybrid import LUT seed: `{not args.hybrid_no_import_lut}`")
        lines.append(f"- Hybrid stage guard: `{not args.hybrid_no_stage_guard}`")
    lines.extend([
        f"- Restarts: {args.restarts}",
        f"- Steps per restart: {args.steps}",
        f"- Bailout steps: {args.bailout_steps}",
        f"- Decompressors: `{args.decompressors}`",
        f"- Compressors: `{args.compressors}`",
        f"- Raw CSV: `{csv_path}`",
        f"- Best CSV: `{best_csv}`",
        f"- Selected artifacts: `{selected_dir}`",
        f"- Total step rows: {len(rows)}",
        f"- Valid best benchmarks: {len(valid_best)}",
        f"- Best improved/same/worse vs baseline `{args.baseline_compressor}`: {improved}/{same}/{worse}",
        f"- Best geomean ratios vs baseline: gates={fmt(gate_geo)}, depth={fmt(depth_geo)}, inv={fmt(inv_geo)}, wall={fmt(wall_geo)}",
        "",
        "| benchmark | best compressor | decompressor | gates | depth | inv | restart | step | equiv |",
        "|---|---|---|---:|---:|---:|---:|---:|---|",
    ])
    for row in sorted(best_rows, key=lambda r: str(r["benchmark"])):
        lines.append(
            f"| {row['benchmark']} | `{row.get('compressor', '')}` | `{row.get('decompressor', '')}` | "
            f"{row.get('metric_gates', '')} | {row.get('metric_depth', '')} | {row.get('metric_inv_edges', '')} | "
            f"{row.get('restart', '')} | {row.get('step', '')} | {row.get('equiv_status', '')} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Randomized MIG/mMIG design-space exploration inspired by Lee et al. DAC 2024 LBR.")
    parser.add_argument("--binary", default="./build_check/examples/blif2mig_2")
    parser.add_argument("--equ-binary", default="./build_check/examples/equ")
    parser.add_argument("--abc-bin", default="abc")
    parser.add_argument("--bench-dir", required=True)
    parser.add_argument("--bench", action="append", default=[])
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--subject", choices=["mmig", "mig", "hybrid"], default="mmig")
    parser.add_argument("--mode", choices=["area", "depth"], default="area")
    parser.add_argument("--mig-flow", default="dac19_compat")
    parser.add_argument("--mmig-max-iters", type=int, default=3)
    parser.add_argument("--mmig-advanced-rounds", type=int, default=2)
    parser.add_argument("--mmig-stage", choices=["pre", "post", "both"], default="both")
    parser.add_argument("--hybrid-mmig-mig-flow", default="identity")
    parser.add_argument("--hybrid-no-import-lut", action="store_true", default=False)
    parser.add_argument("--hybrid-no-stage-guard", action="store_true", default=False)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--bailout-steps", type=int, default=4)
    parser.add_argument("--restart-timeout-s", type=int, default=0)
    parser.add_argument("--step-timeout-s", type=int, default=600)
    parser.add_argument("--abc-timeout-s", type=int, default=120)
    parser.add_argument("--equiv-timeout-s", type=int, default=180)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--decompressors",
        default=(
            "identity,abc_resyn2,abc_resyn2rs,abc_if3,abc_if4,abc_if5,abc_if6,abc_if6_mfs,"
            "abc_dch_if3,abc_dch_if4,abc_dch_if5,abc_dch_if6,abc_dch_if6_mfs,abc_dc2,abc_fraig"
        ),
    )
    parser.add_argument("--compressors", default="")
    parser.add_argument("--baseline-compressor", default="")
    parser.add_argument("--no-baseline", action="store_true", default=False)
    parser.add_argument("--no-equivalence", action="store_true", default=False)
    parser.add_argument("--resume", action="store_true", default=False)
    parser.add_argument(
        "--inv-weight",
        type=float,
        default=0.0,
        help=(
            "Weight applied to inverter count in the primary objective. "
            "0.0 (default) keeps strict lex order (gates, depth, inv). "
            "Values like 0.25 / 0.5 / 1.0 fold inverters into the primary "
            "key as cost = gates + w*inverters, surfacing mMIG-favored "
            "solutions that pure-MIG DSE cannot rank."
        ),
    )
    args = parser.parse_args()
    set_inv_weight(args.inv_weight)
    if not args.compressors:
        args.compressors = default_compressors(args.subject)
    if not args.baseline_compressor:
        args.baseline_compressor = default_baseline_compressor(args.subject)

    if not Path(args.binary).exists():
        print(f"[error] binary not found: {args.binary}", file=sys.stderr)
        return 2
    bench_dir = Path(args.bench_dir)
    if not bench_dir.exists():
        print(f"[error] benchmark directory not found: {bench_dir}", file=sys.stderr)
        return 2
    benches = resolve_benches(bench_dir, args.bench)
    if not benches:
        print("[error] no benchmarks selected", file=sys.stderr)
        return 2

    compressors = compressor_catalog(args.subject, max(1, args.mmig_advanced_rounds))
    if args.baseline_compressor not in compressors:
        print(f"[error] unknown baseline compressor: {args.baseline_compressor}", file=sys.stderr)
        return 2
    for name in parse_list(args.compressors):
        if name not in compressors:
            print(f"[error] unknown compressor: {name}", file=sys.stderr)
            return 2

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    csv_path = out_root / "dse_raw.csv"
    best_csv = out_root / "dse_best.csv"
    summary_path = out_root / "DSE_SUMMARY.md"
    selected_dir = out_root / "selected"

    all_rows: List[Dict[str, object]] = []
    best_rows: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
        futures = {executor.submit(run_benchmark, (bench, args)): bench for bench in benches}
        for future in as_completed(futures):
            bench = futures[future]
            result = future.result()
            rows = result["rows"]
            best = result["best"]
            all_rows.extend(rows)
            best_rows.append(best)
            status = "ok" if valid_candidate(best) else "failed"
            resumed = " resumed" if result.get("resumed") else ""
            print(
                f"[dse] {bench.name:<18} {status:<7} best={best.get('metric_gates', '')}/{best.get('metric_depth', '')} "
                f"{best.get('compressor', '')}/{best.get('decompressor', '')}{resumed}",
                flush=True,
            )

    all_rows.sort(key=lambda r: (str(r["benchmark"]), int(r.get("restart", 0)), int(r.get("step", 0))))
    best_rows.sort(key=lambda r: str(r["benchmark"]))
    write_csv(csv_path, all_rows)
    write_csv(best_csv, best_rows)
    copy_best_artifacts(selected_dir, best_rows)
    write_summary(summary_path, all_rows, best_rows, args, csv_path, best_csv, selected_dir)
    print(f"[done] raw CSV: {csv_path}")
    print(f"[done] best CSV: {best_csv}")
    print(f"[done] summary: {summary_path}")
    print(f"[done] selected artifacts: {selected_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
