#!/usr/bin/env python3
"""Plot the sinusoidal benchmark results with matplotlib.

Output: 3 (p = 2, 3, 4) x 3 (L2, H1, H2) = 9 figures

    program/figures/sinusoidal_plots/sinusoidal_p2_L2.png
    program/figures/sinusoidal_plots/sinusoidal_p2_H1.png
    ...

Each figure overlays all configurations:
  * RKGSI, q = 0, 1, 2
  * all Gauss variants found in program/results/sinusoidal_gauss_*.csv

The style of every figure can be edited independently in FIGURE_OVERRIDES
below.  Global series styles are in SERIES_STYLES / GAUSS_STYLES.

Dependencies:
    python3 -m pip install matplotlib
    # or: uv run --with matplotlib python program/plot_sinusoidal.py
"""

from __future__ import annotations

import csv
import glob
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "program" / "results"
OUT_DIR = ROOT / "program" / "figures" / "sinusoidal_plots"

PS = (2, 3, 4)
NORMS = ("L2", "H1", "H2")

FIG_W = 4.0       # 3 figures across  -> 12 in
FIG_H = 3.6       # 3 figures stacked -> 10.8 in
DPI = 300

YLAB = {
    "L2": r"$\log_{10}\Vert e\Vert_{L_2}$",
    "H1": r"$\log_{10}\Vert e\Vert_{H_1}$",
    "H2": r"$\log_{10}\Vert e\Vert_{H_2}$",
}
NORM_SUB = {"L2": "L_2", "H1": "H_1", "H2": "H_2"}

# ===========================================================================
# Style configuration -- edit this block
# ===========================================================================

# Default style of one figure.  Every field can be overridden per figure in
# FIGURE_OVERRIDES.
DEFAULT_FIG_STYLE = {
    "figsize": (FIG_W, FIG_H),
    "legend_loc": "best",
    "legend_fontsize": 7,
    "handlelength": 3.0,
    "handletextpad": 0.6,
    "borderaxespad": 0.4,
    "labelspacing": 0.35,
    "title": None,          # None -> "$p = ..., H_i error$"
    "xlim": None,
    "ylim": None,
    "theory_style": {"color": "black", "linestyle": "--", "linewidth": 1.0},
}

# Per-figure overrides.  Keys are (p, norm).  Any DEFAULT_FIG_STYLE field can
# be overridden, and "series" overrides one curve at a time.
#
# Example:
#   (2, "L2"): {
#       "legend_loc": "upper right",
#       "series": {
#           "RKGSI_q0": {"color": "red", "linestyle": "--"},
#           "GI_1": {"marker": "P"},
#       },
#   },
FIGURE_OVERRIDES = {
    # (2, "L2"): {},
    # (2, "H1"): {},
    # (2, "H2"): {},
    # (3, "L2"): {},
    # (3, "H1"): {},
    # (3, "H2"): {},
    # (4, "L2"): {},
    # (4, "H1"): {},
    # (4, "H2"): {},
}

# Base styles of the RKGSI curves.
RKGSI_STYLES = {
    "RKGSI_q0": {
        "label": r"$\mathrm{RKGSI},\ q=0$",
        "color": "C0", "linestyle": "-", "marker": "o",
    },
    "RKGSI_q1": {
        "label": r"$\mathrm{RKGSI},\ q=1$",
        "color": "C1", "linestyle": "-", "marker": "s",
    },
    "RKGSI_q2": {
        "label": r"$\mathrm{RKGSI},\ q=2$",
        "color": "C2", "linestyle": "-", "marker": "^",
    },
}

# Base styles of the Gauss variants.  GI_1, GI_2, ... are assigned in
# chronological order of the matching CSV files.
GAUSS_STYLES = [
    {"label": r"$\mathrm{GI}_{1}$", "color": "C3", "linestyle": "--", "marker": "D"},
    {"label": r"$\mathrm{GI}_{2}$", "color": "C4", "linestyle": ":",  "marker": "x"},
    {"label": r"$\mathrm{GI}_{3}$", "color": "C5", "linestyle": "-.", "marker": "v"},
    {"label": r"$\mathrm{GI}_{4}$", "color": "C6", "linestyle": (0, (3, 1, 1, 1)), "marker": "P"},
    {"label": r"$\mathrm{GI}_{5}$", "color": "C7", "linestyle": (0, (1, 1)), "marker": "*"},
]

# ===========================================================================


def _float(row: dict, key: str):
    value = row.get(key, "")
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_results(path: Path):
    rows = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("valid", "yes") not in ("yes", ""):
                continue
            p = _float(row, "p")
            h = _float(row, "h")
            if p is None or h is None or h <= 0.0:
                continue
            vals = {name: _float(row, name) for name in NORMS}
            if any(v is None or v <= 0.0 for v in vals.values()):
                continue
            q = _float(row, "q")
            rows.append({
                "p": int(p),
                "q": None if q is None else int(q),
                "h": h,
                **vals,
            })
    return rows


def latest_rkgsi():
    candidates = sorted(
        glob.glob(str(RESULT_DIR / "sinusoidal_rkgsi_*.csv")),
        key=os.path.getmtime,
    )
    candidates = [p for p in candidates if not p.endswith(("_rates.csv", "_config.csv"))]
    fallback = RESULT_DIR / "sinusoidal_results.csv"
    if candidates:
        return Path(candidates[-1])
    if fallback.is_file():
        return fallback
    return None


def gauss_runs():
    files = sorted(
        [p for p in glob.glob(str(RESULT_DIR / "sinusoidal_gauss_*.csv"))
         if not p.endswith(("_rates.csv", "_config.csv"))],
        key=os.path.getmtime,
    )
    runs = []
    for i, path in enumerate(files):
        base = GAUSS_STYLES[i] if i < len(GAUSS_STYLES) else {
            "label": rf"$\mathrm{{GI}}_{{{i + 1}}}$",
            "color": f"C{i % 10}", "linestyle": "--", "marker": "o",
        }
        runs.append((f"GI_{i + 1}", Path(path), dict(base)))
    return runs


def merged_fig_style(p: int, norm: str):
    style = dict(DEFAULT_FIG_STYLE)
    style.update(FIGURE_OVERRIDES.get((p, norm), {}))
    return style


def merged_series_style(base: dict, fig_style: dict, series_id: str):
    style = dict(base)
    overrides = fig_style.get("series", {}).get(series_id, {})
    style.update(overrides)
    return style


def theory_slope(p: int, norm: str) -> int:
    """Optimal RKGSI rate for the given p and error norm."""
    if norm == "L2":
        return 2 if p == 2 else p + 1
    if norm == "H1":
        return p
    return p - 1


def add_reference_line(ax, slope: int, label: str, style: dict):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xa = x0 + 0.05*(x1 - x0)
    xb = x0 + 0.30*(x1 - x0)
    yb = y1 - 0.10*(y1 - y0)
    ya = yb - slope*(xb - xa)
    if ya < y0 + 0.05*(y1 - y0):
        ya = y0 + 0.05*(y1 - y0)
        yb = ya + slope*(xb - xa)
    ax.plot([xa, xb], [ya, yb], label=label, zorder=1, **style)


def plot_one(p: int, norm: str, rkgsi, gauss):
    fs = merged_fig_style(p, norm)
    fig, ax = plt.subplots(figsize=fs["figsize"])

    # RKGSI q = 0, 1, 2
    for q in (0, 1, 2):
        pts = sorted(
            [r for r in rkgsi if r["p"] == p and r["q"] == q],
            key=lambda r: r["h"],
        )
        if not pts:
            continue
        sid = f"RKGSI_q{q}"
        st = merged_series_style(RKGSI_STYLES[sid], fs, sid)
        x = [math.log10(r["h"]) for r in pts]
        y = [math.log10(r[norm]) for r in pts]
        ax.plot(x, y, label=st["label"], color=st["color"],
                linestyle=st["linestyle"], marker=st["marker"],
                markersize=4.5, linewidth=1.3)

    # Gauss variants
    for sid, path, base in gauss:
        pts = sorted(
            [r for r in read_results(path) if r["p"] == p],
            key=lambda r: r["h"],
        )
        if not pts:
            continue
        st = merged_series_style(base, fs, sid)
        x = [math.log10(r["h"]) for r in pts]
        y = [math.log10(r[norm]) for r in pts]
        ax.plot(x, y, label=st["label"], color=st["color"],
                linestyle=st["linestyle"], marker=st["marker"],
                markersize=4.5, linewidth=1.3)

    ax.set_xlabel(r"$\log_{10} h$")
    ax.set_ylabel(YLAB[norm])
    title = fs["title"] or rf"$p = {p},\ {NORM_SUB[norm]}\ \mathrm{{error}}$"
    ax.set_title(title)
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.7)
    add_reference_line(ax, theory_slope(p, norm),
                       rf"$\mathrm{{theory}}: {theory_slope(p, norm)}$",
                       fs["theory_style"])
    if fs["xlim"] is not None:
        ax.set_xlim(fs["xlim"])
    if fs["ylim"] is not None:
        ax.set_ylim(fs["ylim"])
    ax.legend(fontsize=fs["legend_fontsize"], frameon=False,
              loc=fs["legend_loc"], handlelength=fs["handlelength"],
              handletextpad=fs["handletextpad"],
              borderaxespad=fs["borderaxespad"],
              labelspacing=fs["labelspacing"])
    fig.tight_layout()

    out = OUT_DIR / f"sinusoidal_p{p}_{norm}.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rkgsi_path = latest_rkgsi()
    if rkgsi_path is None:
        raise SystemExit("no RKGSI CSV found in program/results/")
    rkgsi = read_results(rkgsi_path)
    gauss = gauss_runs()

    print(f"RKGSI: {rkgsi_path}")
    for sid, path, _ in gauss:
        print(f"Gauss: {path}  ->  {sid}")

    written = []
    for p in PS:
        for norm in NORMS:
            written.append(plot_one(p, norm, rkgsi, gauss))

    print("written:")
    for path in written:
        print(" ", path)


if __name__ == "__main__":
    main()
