#!/usr/bin/env python3
"""Generate fig/mmig_comparison.pdf and fig/mmig_comparison.png from data/results.json.

Run from the repo root:
    python3 scripts/gen_figure.py
"""

import json
import os
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "1700000000")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

repo = Path(__file__).resolve().parents[1]
data = json.loads((repo / "data" / "results.json").read_text())

rows = sorted(data["rows"], key=lambda r: r[13])
benchmarks = [r[0] for r in rows]
dG = [-r[12] for r in rows]
dI = [-r[13] for r in rows]
dA = [-r[14] for r in rows]

y = np.arange(len(benchmarks))

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(3.5, 3.4), sharey=True,
    gridspec_kw={"width_ratios": [1.55, 1], "wspace": 0.06},
)

h, h2 = 0.50, 0.22

colors_dI = ["#2E7D32" if v > 0 else "#BDBDBD" for v in dI]
ax1.barh(y, dI, h, color=colors_dI, zorder=3)
ax1.axvline(0, color="black", linewidth=0.6, zorder=2)
ax1.set_xlim(-1, 39)
ax1.set_xlabel(r"$\Delta I$ (%)", fontsize=7.5)
ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
ax1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}"))
ax1.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.6, zorder=1)
ax1.set_yticks(y)
ax1.set_yticklabels(benchmarks, fontsize=7.5)
ax1.tick_params(axis="both", labelsize=7.5)
ax1.invert_yaxis()

for i, v in enumerate(dI):
    if v > 2:
        ax1.text(v + 0.6, i, f"{v:.0f}%", va="center", fontsize=6.5, color="#1B5E20")

ax1.text(0.97, 0.02, r"$\Delta I$ (inverters)", transform=ax1.transAxes,
         ha="right", va="bottom", fontsize=6.5,
         bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))

ax2.barh(y - 0.13, dG, h2, color="#1565C0", zorder=3, label=r"$\Delta G$")
ax2.barh(y + 0.13, dA, h2, color="#C62828", zorder=3, label=r"$\Delta A$")
ax2.axvline(0, color="black", linewidth=0.6, zorder=2)
ax2.set_xlim(-0.3, 10)
ax2.set_xlabel(r"$\Delta G$, $\Delta A$ (%)", fontsize=7.5)
ax2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
ax2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}"))
ax2.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.6, zorder=1)
ax2.tick_params(axis="both", labelsize=7.5)
ax2.legend(fontsize=6.5, loc="lower right", ncol=2, framealpha=0.85,
           handlelength=1.0, handletextpad=0.3, columnspacing=0.5, borderpad=0.3)

fig.subplots_adjust(left=0.26, right=0.99, top=0.99, bottom=0.11, wspace=0.06)

out_dir = repo / "fig"
out_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(
    out_dir / "mmig_comparison.pdf",
    bbox_inches="tight",
    dpi=200,
    metadata={"Creator": "scripts/gen_figure.py", "Producer": "matplotlib"},
)
fig.savefig(out_dir / "mmig_comparison.png", bbox_inches="tight", dpi=200)
print("Saved fig/mmig_comparison.pdf and fig/mmig_comparison.png")
