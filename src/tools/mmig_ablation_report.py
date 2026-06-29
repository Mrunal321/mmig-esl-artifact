#!/usr/bin/env python3

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from epfl_flow_report import run_case


@dataclass(frozen=True)
class Profile:
    name: str
    description: str
    use_mmig: bool
    mmig_flow: str = "round_robin"
    extra_args: Tuple[str, ...] = ()


def profile_catalog(rounds: int) -> Dict[str, Profile]:
    adv = ("--mmig-advanced", f"--mmig-advanced-rounds={rounds}")
    paper_adv = ("--mmig-advanced", "--mmig-advanced-rounds=1")
    legacy = (
        "--no-mmig-rewrite-ranking",
        "--no-mmig-normalize-inner",
        "--no-mmig-resub-ranking",
        "--no-mmig-auto-sr5",
        "--mmig-adv-no-tuned-policy",
    )
    return {
        "mig_dac19": Profile(
            "mig_dac19",
            "MIG-only dac19_default baseline; no mMIG stage",
            False,
        ),
        "rr_legacy": Profile(
            "rr_legacy",
            "round_robin with new ranking/normalization/resub ranking/auto-SR5/tuned-policy disabled",
            True,
            "round_robin",
            adv + legacy,
        ),
        "rr_rank": Profile(
            "rr_rank",
            "round_robin with candidate ranking enabled, complemented-inner normalization disabled",
            True,
            "round_robin",
            adv + ("--no-mmig-normalize-inner", "--no-mmig-auto-sr5", "--mmig-adv-no-tuned-policy"),
        ),
        "rr_norm": Profile(
            "rr_norm",
            "round_robin with complemented-inner normalization enabled, candidate ranking disabled",
            True,
            "round_robin",
            adv + ("--no-mmig-rewrite-ranking", "--no-mmig-resub-ranking", "--no-mmig-auto-sr5", "--mmig-adv-no-tuned-policy"),
        ),
        "rr_rank_norm": Profile(
            "rr_rank_norm",
            "round_robin with current rewrite ranking and complemented-inner normalization",
            True,
            "round_robin",
            adv,
        ),
        "compress2rs": Profile(
            "compress2rs",
            "advanced compress2rs profile without dont-cares or interleaved seeding",
            True,
            "compress2rs",
            adv,
        ),
        "compress2rs_dc": Profile(
            "compress2rs_dc",
            "advanced compress2rs with local don't-cares",
            True,
            "compress2rs",
            adv + ("--mmig-dont-cares",),
        ),
        "compress2rs_dc_reseed": Profile(
            "compress2rs_dc_reseed",
            "advanced compress2rs with local don't-cares and interleaved minority seeding",
            True,
            "compress2rs",
            adv + ("--mmig-dont-cares", "--mmig-interleaved-seeding"),
        ),
        "paper2019": Profile(
            "paper2019",
            "exact Section 3.1 generic resynthesis sequence from the 2019 paper",
            True,
            "paper2019",
            paper_adv + ("--mmig-adv-no-exact", "--mmig-adv-no-balance", "--mmig-adv-no-objective-guard", "--mmig-adv-no-stagnation-stop"),
        ),
        "paper2019_guarded": Profile(
            "paper2019_guarded",
            "paper2019 sequence with the mMIG round objective guard left enabled",
            True,
            "paper2019",
            paper_adv + ("--mmig-adv-no-exact", "--mmig-adv-no-balance"),
        ),
        "epfl_dc_reseed": Profile(
            "epfl_dc_reseed",
            "heavier EPFL profile with local don't-cares and interleaved minority seeding",
            True,
            "epfl",
            adv + ("--mmig-dont-cares", "--mmig-interleaved-seeding"),
        ),
    }


PROFILE_SETS = {
    "publication_seed": (
        "rr_legacy",
        "rr_rank",
        "rr_norm",
        "rr_rank_norm",
        "compress2rs",
        "compress2rs_dc",
        "compress2rs_dc_reseed",
    ),
    "portfolio": ("rr_rank_norm", "compress2rs_dc_reseed"),
    "quick": ("rr_legacy", "rr_rank_norm", "compress2rs_dc_reseed"),
    "paper_flow": ("rr_legacy", "paper2019", "paper2019_guarded", "compress2rs_dc_reseed"),
    "full": (
        "mig_dac19",
        "rr_legacy",
        "rr_rank",
        "rr_norm",
        "rr_rank_norm",
        "compress2rs",
        "compress2rs_dc",
        "compress2rs_dc_reseed",
        "paper2019",
        "paper2019_guarded",
        "epfl_dc_reseed",
    ),
}


def resolve_benches(bench_dir: Path, selectors: Sequence[str]) -> List[Path]:
    benches = sorted(bench_dir.glob("*.blif"))
    if not selectors:
        return benches

    available = {b.name: b for b in benches}
    selected = []
    seen = set()
    for item in selectors:
        base = Path(item).name
        names = [base] if base.endswith(".blif") else [f"{base}.blif", f"{base}_rebuilt.blif"]
        match = next((name for name in names if name in available), None)
        if match is None:
            print(f"[warn] benchmark selector did not match: {item}", file=sys.stderr)
            continue
        if match not in seen:
            selected.append(available[match])
            seen.add(match)
    return selected


def parse_profiles(value: str, catalog: Dict[str, Profile]) -> List[Profile]:
    names: List[str] = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        if token in PROFILE_SETS:
            names.extend(PROFILE_SETS[token])
        else:
            names.append(token)

    profiles = []
    seen = set()
    for name in names:
        if name not in catalog:
            valid = ", ".join(sorted(catalog))
            raise ValueError(f"unknown profile '{name}'. Valid profiles: {valid}")
        if name not in seen:
            profiles.append(catalog[name])
            seen.add(name)
    return profiles


def as_int(row: Dict[str, object], key: str, missing: int = 10**12) -> int:
    value = row.get(key, "")
    if value in ("", None):
        return missing
    return int(value)


def as_float(row: Dict[str, object], key: str) -> Optional[float]:
    value = row.get(key, "")
    if value in ("", None):
        return None
    return float(value)


def objective(row: Dict[str, object], mode: str) -> Tuple[int, int, int]:
    if mode == "area":
        return (as_int(row, "metric_gates"), as_int(row, "metric_depth"), as_int(row, "metric_inv_edges"))
    return (as_int(row, "metric_depth"), as_int(row, "metric_gates"), as_int(row, "metric_inv_edges"))


def geometric_mean(values: Iterable[float]) -> Optional[float]:
    positive = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positive:
        return None
    return math.exp(sum(math.log(v) for v in positive) / len(positive))


def run_equivalence(equ_binary: Path, bench: Path, revised: Path, timeout_s: int) -> Tuple[str, str]:
    if not equ_binary.exists():
        return "skipped", f"equivalence binary not found: {equ_binary}"
    if not revised.exists():
        return "missing_output", f"optimized BLIF not found: {revised}"
    try:
        completed = subprocess.run(
            [str(equ_binary), str(bench), str(revised)],
            capture_output=True,
            text=True,
            timeout=timeout_s if timeout_s > 0 else None,
        )
    except subprocess.TimeoutExpired:
        return "timeout", f"equivalence check exceeded {timeout_s} seconds"

    output = ((completed.stdout or "") + ("\n" + completed.stderr if completed.stderr else "")).strip()
    if completed.returncode == 0 and "EQUIVALENT" in output:
        return "equivalent", output
    if completed.returncode == 2:
        return "not_equivalent", output
    return f"error_{completed.returncode}", output


def make_row(
    profile: Profile,
    bench: Path,
    mode: str,
    mig_flow: str,
    result: Dict[str, object],
    log_path: Path,
    out_base: Path,
    equiv_status: str,
    equiv_output: str,
) -> Dict[str, object]:
    mig = result["mig_metrics"]
    mmig = result["mmig_metrics"]
    cec = result["cec"]
    metric_source = "mmig" if profile.use_mmig else "mig"
    metric = mmig if profile.use_mmig else mig
    return {
        "benchmark": bench.name,
        "profile": profile.name,
        "description": profile.description,
        "mode": mode,
        "mig_flow": mig_flow,
        "mmig_flow": profile.mmig_flow if profile.use_mmig else "",
        "metric_source": metric_source,
        "return_code": result["return_code"],
        "wall_s": f"{float(result['wall_s']):.3f}",
        "mig_time_s": result["mig_time_s"],
        "mig_gates": mig["gates"],
        "mig_depth": mig["depth"],
        "mig_majority_nodes": mig["majority_nodes"],
        "mig_minority_nodes": mig["minority_nodes"],
        "mig_inv_edges": mig["inv_edges_total"],
        "mmig_time_s": result["mmig_time_s"] if profile.use_mmig else "",
        "mmig_gates": mmig["gates"] if profile.use_mmig else "",
        "mmig_depth": mmig["depth"] if profile.use_mmig else "",
        "mmig_majority_nodes": mmig["majority_nodes"] if profile.use_mmig else "",
        "mmig_minority_nodes": mmig["minority_nodes"] if profile.use_mmig else "",
        "mmig_inv_edges": mmig["inv_edges_total"] if profile.use_mmig else "",
        "metric_gates": metric["gates"],
        "metric_depth": metric["depth"],
        "metric_inv_edges": metric["inv_edges_total"],
        "cec_accepted": cec["cec_accepted"] if profile.use_mmig else "",
        "cec_rejected": cec["cec_rejected"] if profile.use_mmig else "",
        "cec_inconclusive": cec["cec_inconclusive"] if profile.use_mmig else "",
        "cec_time_ms": cec["cec_time_ms"] if profile.use_mmig else "",
        "stage_reject_lines": result["stage_reject_lines"] if profile.use_mmig else "",
        "equiv_status": equiv_status,
        "equiv_output": equiv_output,
        "log_path": str(log_path),
        "out_base": str(out_base),
        "command": result["command"],
    }


def run_task(task: Dict[str, object]) -> Dict[str, object]:
    profile: Profile = task["profile"]  # type: ignore[assignment]
    bench: Path = task["bench"]  # type: ignore[assignment]
    out_root: Path = task["out_root"]  # type: ignore[assignment]
    mode: str = task["mode"]  # type: ignore[assignment]
    mig_flow: str = task["mig_flow"]  # type: ignore[assignment]
    binary: Path = task["binary"]  # type: ignore[assignment]
    equ_binary: Path = task["equ_binary"]  # type: ignore[assignment]
    timeout_s: int = task["timeout_s"]  # type: ignore[assignment]
    equiv_timeout_s: int = task["equiv_timeout_s"]  # type: ignore[assignment]
    mmig_max_iters: int = task["mmig_max_iters"]  # type: ignore[assignment]
    mmig_stage: str = task["mmig_stage"]  # type: ignore[assignment]
    run_equiv: bool = task["run_equiv"]  # type: ignore[assignment]
    resume: bool = task["resume"]  # type: ignore[assignment]

    out_base = out_root / f"{bench.stem}_{profile.name}"
    log_path = out_root / f"{bench.stem}_{profile.name}.log"
    status_path = out_root / f"{bench.stem}_{profile.name}.json"
    if resume and status_path.exists() and log_path.exists():
        with status_path.open(encoding="utf-8") as fp:
            row = json.load(fp)
        row["resumed"] = "1"
        return row

    result = run_case(
        binary=binary,
        bench=bench,
        out_base=out_base,
        mode=mode,
        mig_flow=mig_flow,
        use_mmig=profile.use_mmig,
        mmig_flow=profile.mmig_flow,
        mmig_max_iters=mmig_max_iters,
        mmig_stage=mmig_stage,
        mmig_cec=True,
        use_zero_gain=False,
        extra_args=list(profile.extra_args),
        timeout_s=timeout_s,
    )
    log_path.write_text(result["stdout"], encoding="utf-8", errors="ignore")

    equiv_status = "skipped"
    equiv_output = ""
    if run_equiv and int(result["return_code"]) == 0:
        revised = Path(str(out_base) + ("_mmig_opt.blif" if profile.use_mmig else "_maj_opt.blif"))
        equiv_status, equiv_output = run_equivalence(equ_binary, bench, revised, equiv_timeout_s)

    row = make_row(profile, bench, mode, mig_flow, result, log_path, out_base, equiv_status, equiv_output)
    row["resumed"] = "0"
    with status_path.open("w", encoding="utf-8") as fp:
        json.dump(row, fp, indent=2, sort_keys=True)
    return row


FIELDNAMES = [
    "benchmark",
    "profile",
    "description",
    "mode",
    "mig_flow",
    "mmig_flow",
    "metric_source",
    "return_code",
    "wall_s",
    "mig_time_s",
    "mig_gates",
    "mig_depth",
    "mig_majority_nodes",
    "mig_minority_nodes",
    "mig_inv_edges",
    "mmig_time_s",
    "mmig_gates",
    "mmig_depth",
    "mmig_majority_nodes",
    "mmig_minority_nodes",
    "mmig_inv_edges",
    "metric_gates",
    "metric_depth",
    "metric_inv_edges",
    "cec_accepted",
    "cec_rejected",
    "cec_inconclusive",
    "cec_time_ms",
    "stage_reject_lines",
    "equiv_status",
    "equiv_output",
    "log_path",
    "out_base",
    "command",
    "resumed",
]


def successful(row: Dict[str, object]) -> bool:
    return (
        str(row.get("return_code")) == "0"
        and row.get("metric_gates") not in ("", None)
        and row.get("metric_depth") not in ("", None)
        and row.get("metric_inv_edges") not in ("", None)
        and row.get("equiv_status") in ("equivalent", "skipped")
    )


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def copy_selected_artifacts(selected_dir: Path, selected: List[Dict[str, object]]) -> None:
    selected_dir.mkdir(parents=True, exist_ok=True)
    for row in selected:
        out_base = Path(str(row["out_base"]))
        dst_base = selected_dir / f"{Path(str(row['benchmark'])).stem}_{row['profile']}"
        for suffix in ("_mmig_opt.blif", "_mmig_opt.v", "_maj_opt.blif", "_maj_opt.v"):
            src = Path(str(out_base) + suffix)
            if src.exists():
                shutil.copy2(src, Path(str(dst_base) + suffix))
        log_path = Path(str(row["log_path"]))
        if log_path.exists():
            shutil.copy2(log_path, Path(str(dst_base) + ".log"))


def compare_rows(candidate: Dict[str, object], baseline: Dict[str, object], mode: str) -> int:
    cand_obj = objective(candidate, mode)
    base_obj = objective(baseline, mode)
    if cand_obj < base_obj:
        return -1
    if cand_obj > base_obj:
        return 1
    return 0


def ratio(candidate: Dict[str, object], baseline: Dict[str, object], key: str) -> Optional[float]:
    c = as_float(candidate, key)
    b = as_float(baseline, key)
    if c is None or b is None or b <= 0.0:
        return None
    return c / b


def format_ratio(value: Optional[float]) -> str:
    if value is None:
        return "n/a"
    return f"{value:.4f}"


def write_summary(
    path: Path,
    rows: List[Dict[str, object]],
    profiles: List[Profile],
    baseline_profile: str,
    mode: str,
    csv_path: Path,
    selected_csv: Path,
    selected_dir: Path,
) -> None:
    by_bench: Dict[str, Dict[str, Dict[str, object]]] = {}
    for row in rows:
        by_bench.setdefault(str(row["benchmark"]), {})[str(row["profile"])] = row

    selected = []
    for bench, variants in sorted(by_bench.items()):
        candidates = [row for row in variants.values() if successful(row)]
        if not candidates:
            continue
        selected.append(min(candidates, key=lambda row: objective(row, mode)))

    baseline_pairs = []
    improved = same = worse = 0
    for row in selected:
        baseline = by_bench[str(row["benchmark"])].get(baseline_profile)
        if baseline is None or not successful(baseline):
            continue
        baseline_pairs.append((row, baseline))
        cmp_value = compare_rows(row, baseline, mode)
        if cmp_value < 0:
            improved += 1
        elif cmp_value == 0:
            same += 1
        else:
            worse += 1

    selected_geo = {
        "gates": geometric_mean(r for r in (ratio(row, base, "metric_gates") for row, base in baseline_pairs) if r is not None),
        "depth": geometric_mean(r for r in (ratio(row, base, "metric_depth") for row, base in baseline_pairs) if r is not None),
        "inv": geometric_mean(r for r in (ratio(row, base, "metric_inv_edges") for row, base in baseline_pairs) if r is not None),
        "wall": geometric_mean(r for r in (ratio(row, base, "wall_s") for row, base in baseline_pairs) if r is not None),
    }

    lines = [
        "# mMIG EPFL Ablation Summary",
        "",
        f"- Mode: `{mode}`",
        f"- Baseline profile: `{baseline_profile}`",
        f"- Raw CSV: `{csv_path}`",
        f"- Selected CSV: `{selected_csv}`",
        f"- Selected artifacts: `{selected_dir}`",
        f"- Runs: {len(rows)}",
        f"- Successful/equivalent runs: {sum(1 for row in rows if successful(row))}",
        f"- Selected benchmarks: {len(selected)}",
        f"- Selected improved/same/worse vs `{baseline_profile}` by objective: {improved}/{same}/{worse}",
        f"- Selected geomean ratios vs `{baseline_profile}`: gates={format_ratio(selected_geo['gates'])}, depth={format_ratio(selected_geo['depth'])}, inv={format_ratio(selected_geo['inv'])}, wall={format_ratio(selected_geo['wall'])}",
        "",
        "## Profiles",
        "",
        "| profile | description |",
        "|---|---|",
    ]
    for profile in profiles:
        lines.append(f"| `{profile.name}` | {profile.description} |")

    lines.extend(
        [
            "",
            "## Profile Comparison Vs Baseline",
            "",
            "| profile | common | improved | same | worse | gmean gates | gmean depth | gmean inv | gmean wall | wins |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )

    win_counts = {profile.name: 0 for profile in profiles}
    for row in selected:
        win_counts[str(row["profile"])] = win_counts.get(str(row["profile"]), 0) + 1

    for profile in profiles:
        pairs = []
        p_improved = p_same = p_worse = 0
        for bench, variants in by_bench.items():
            row = variants.get(profile.name)
            baseline = variants.get(baseline_profile)
            if row is None or baseline is None or not successful(row) or not successful(baseline):
                continue
            pairs.append((row, baseline))
            cmp_value = compare_rows(row, baseline, mode)
            if cmp_value < 0:
                p_improved += 1
            elif cmp_value == 0:
                p_same += 1
            else:
                p_worse += 1
        geos = {
            "gates": geometric_mean(r for r in (ratio(row, base, "metric_gates") for row, base in pairs) if r is not None),
            "depth": geometric_mean(r for r in (ratio(row, base, "metric_depth") for row, base in pairs) if r is not None),
            "inv": geometric_mean(r for r in (ratio(row, base, "metric_inv_edges") for row, base in pairs) if r is not None),
            "wall": geometric_mean(r for r in (ratio(row, base, "wall_s") for row, base in pairs) if r is not None),
        }
        lines.append(
            f"| `{profile.name}` | {len(pairs)} | {p_improved} | {p_same} | {p_worse} | "
            f"{format_ratio(geos['gates'])} | {format_ratio(geos['depth'])} | {format_ratio(geos['inv'])} | {format_ratio(geos['wall'])} | {win_counts.get(profile.name, 0)} |"
        )

    lines.extend(
        [
            "",
            "## Selected Best Per Benchmark",
            "",
            "| benchmark | selected | gates | depth | inv | dg | dd | di | equiv |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in selected:
        baseline = by_bench[str(row["benchmark"])].get(baseline_profile, row)
        dg = as_int(row, "metric_gates") - as_int(baseline, "metric_gates")
        dd = as_int(row, "metric_depth") - as_int(baseline, "metric_depth")
        di = as_int(row, "metric_inv_edges") - as_int(baseline, "metric_inv_edges")
        lines.append(
            f"| {row['benchmark']} | `{row['profile']}` | {row['metric_gates']} | {row['metric_depth']} | {row['metric_inv_edges']} | "
            f"{dg} | {dd} | {di} | {row['equiv_status']} |"
        )

    failures = [row for row in rows if str(row.get("return_code")) != "0" or row.get("equiv_status") not in ("equivalent", "skipped")]
    if failures:
        lines.extend(["", "## Failures Or Non-equivalent Runs", "", "| benchmark | profile | rc | equiv |", "|---|---|---:|---|"])
        for row in failures:
            lines.append(f"| {row['benchmark']} | `{row['profile']}` | {row['return_code']} | {row['equiv_status']} |")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    selected_rows = []
    for row in selected:
        selected_rows.append(dict(row))
    write_csv(selected_csv, selected_rows)
    copy_selected_artifacts(selected_dir, selected_rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run mMIG EPFL profile ablation with resume, parallelism, equivalence, and geomean summaries.")
    parser.add_argument("--binary", default="./build_check/examples/blif2mig_2")
    parser.add_argument("--equ-binary", default="./build_check/examples/equ")
    parser.add_argument("--bench-dir", required=True)
    parser.add_argument("--bench", action="append", default=[])
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--csv", default="")
    parser.add_argument("--selected-csv", default="")
    parser.add_argument("--summary", default="")
    parser.add_argument("--mode", choices=["depth", "area"], default="depth")
    parser.add_argument("--mig-flow", default="dac19_default", choices=["identity", "standard", "dac19_default", "dac19_area", "dac19_compat", "compress2rs", "legacy"])
    parser.add_argument("--profiles", default="publication_seed", help="Comma-separated profile names or sets: quick, portfolio, paper_flow, publication_seed, full")
    parser.add_argument("--baseline-profile", default="rr_legacy")
    parser.add_argument("--mmig-max-iters", type=int, default=3)
    parser.add_argument("--mmig-advanced-rounds", type=int, default=2)
    parser.add_argument("--mmig-stage", choices=["pre", "post", "both"], default="both")
    parser.add_argument("--timeout-s", type=int, default=900)
    parser.add_argument("--equiv-timeout-s", type=int, default=180)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--resume", action="store_true", default=False)
    parser.add_argument("--no-equivalence", action="store_true", default=False)
    args = parser.parse_args()

    binary = Path(args.binary)
    if not binary.exists():
        print(f"[error] binary not found: {binary}", file=sys.stderr)
        return 2
    bench_dir = Path(args.bench_dir)
    if not bench_dir.exists():
        print(f"[error] benchmark directory not found: {bench_dir}", file=sys.stderr)
        return 2

    catalog = profile_catalog(max(1, int(args.mmig_advanced_rounds)))
    try:
        profiles = parse_profiles(args.profiles, catalog)
    except ValueError as err:
        print(f"[error] {err}", file=sys.stderr)
        return 2
    if args.baseline_profile not in [profile.name for profile in profiles]:
        print(f"[error] baseline profile '{args.baseline_profile}' is not in selected profiles", file=sys.stderr)
        return 2

    benches = resolve_benches(bench_dir, args.bench)
    if not benches:
        print("[error] no benchmarks selected", file=sys.stderr)
        return 2

    out_root = Path(args.out_root)
    raw_root = out_root / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)
    csv_path = Path(args.csv) if args.csv else out_root / "ablation_raw.csv"
    selected_csv = Path(args.selected_csv) if args.selected_csv else out_root / "ablation_selected_best.csv"
    summary = Path(args.summary) if args.summary else out_root / "ABLATION_SUMMARY.md"
    selected_dir = out_root / "selected"

    tasks = []
    for bench in benches:
        for profile in profiles:
            tasks.append(
                {
                    "profile": profile,
                    "bench": bench,
                    "out_root": raw_root,
                    "mode": args.mode,
                    "mig_flow": args.mig_flow,
                    "binary": binary,
                    "equ_binary": Path(args.equ_binary),
                    "timeout_s": max(0, int(args.timeout_s)),
                    "equiv_timeout_s": max(0, int(args.equiv_timeout_s)),
                    "mmig_max_iters": max(1, int(args.mmig_max_iters)),
                    "mmig_stage": args.mmig_stage,
                    "run_equiv": not args.no_equivalence,
                    "resume": args.resume,
                }
            )

    rows: List[Dict[str, object]] = []
    jobs = max(1, int(args.jobs))
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_task = {executor.submit(run_task, task): task for task in tasks}
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            profile: Profile = task["profile"]  # type: ignore[assignment]
            bench: Path = task["bench"]  # type: ignore[assignment]
            try:
                row = future.result()
            except Exception as err:
                print(f"[run] {bench.name:<22} {profile.name:<22} failed: {err}", flush=True)
                raise
            rows.append(row)
            status = "ok" if successful(row) else f"rc={row.get('return_code')} equiv={row.get('equiv_status')}"
            resumed = " resumed" if row.get("resumed") == "1" else ""
            print(f"[run] {bench.name:<22} {profile.name:<22} {status}{resumed}", flush=True)

    rows.sort(key=lambda row: (str(row["benchmark"]), str(row["profile"])))
    write_csv(csv_path, rows)
    write_summary(summary, rows, profiles, args.baseline_profile, args.mode, csv_path, selected_csv, selected_dir)

    print(f"[done] raw CSV: {csv_path}")
    print(f"[done] selected CSV: {selected_csv}")
    print(f"[done] summary: {summary}")
    print(f"[done] selected artifacts: {selected_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
