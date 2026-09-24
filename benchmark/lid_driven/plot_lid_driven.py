#!/usr/bin/env python3
"""Plot lid-driven cavity solver profiles against Ghia et al. (1982).

Reads the per-Re centerline CSVs (u_mid_Re<Re>.csv / v_mid_Re<Re>.csv)
produced by extract_lid_driven.py, overlays them on the Ghia-u.csv / Ghia-v.csv
reference tables, and produces two comparison plots (u along y, v along x),
each shown on screen and saved as a 300 dpi PNG.

Usage:
    python3 benchmark/lid_driven/extract_lid_driven.py   # first, to generate the CSVs
    python3 benchmark/lid_driven/plot_lid_driven.py
    python3 benchmark/lid_driven/plot_lid_driven.py --no-show
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

REYNOLDS_NUMBERS = [100, 400, 1000, 3200]

# Categorical palette, fixed order -- validated as CVD-safe for adjacent pairs
# (dataviz skill, references/palette.md), one color per Re so it stays
# consistent between the u-mid and v-mid plots.
COLORS = {
    100: "#2a78d6",   # blue
    400: "#eb6834",   # orange
    1000: "#1baf7a",  # aqua
    3200: "#eda100",  # yellow
}

GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"
MUTED = "#898781"
INK = "#0b0b0b"


def sim_csv_paths(re):
    return SCRIPT_DIR / f"u_mid_Re{re}.csv", SCRIPT_DIR / f"v_mid_Re{re}.csv"


def load_sim_profile(csv_path, pos_col, val_col):
    if not csv_path.exists():
        raise FileNotFoundError(
            f"missing {csv_path.name} -- run extract_lid_driven.py first"
        )
    df = pd.read_csv(csv_path).sort_values(pos_col)
    return df[pos_col].to_numpy(), df[val_col].to_numpy()


def load_ghia(name):
    return pd.read_csv(SCRIPT_DIR / name)


def style_axes(ax):
    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    # for side in ("top", "right"):
        # ax.spines[side].set_visible(False)
    # for side in ("left", "bottom"):
        # ax.spines[side].set_color(AXIS)

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1)

    ax.set_box_aspect(1)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def plot_u_mid(out_png, show):
    ghia = load_ghia("Ghia-u.csv")

    fig, ax = plt.subplots(figsize=(6, 6))
    for re in REYNOLDS_NUMBERS:
        color = COLORS[re]
        u_csv, _ = sim_csv_paths(re)
        y, u = load_sim_profile(u_csv, "Points:1", "vel:0")
        ax.plot(u, y, color=color, linewidth=2, label=f"Re = {re}", zorder=3)
        ax.plot(
            ghia[f"Re={re}"], ghia["y"],
            linestyle="none", marker="o", markersize=6,
            markerfacecolor="none", markeredgecolor=color, markeredgewidth=1.3,
            zorder=4,
        )

    ax.set_xlabel("u")
    ax.set_ylabel("y")
    ax.set_xlim(-0.6, 1.2)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title("u along vertical centerline (x = 0.5)")
    style_axes(ax)
    ax.legend(loc="best", frameon=False)
    fig.text(0.5, 0.01, "line: present solver   o: Ghia et al. (1982)",
              ha="center", color=INK, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")
    if show:
        plt.show(block=False)


def plot_v_mid(out_png, show):
    ghia = load_ghia("Ghia-v.csv")

    fig, ax = plt.subplots(figsize=(6, 6))
    for re in REYNOLDS_NUMBERS:
        color = COLORS[re]
        _, v_csv = sim_csv_paths(re)
        x, v = load_sim_profile(v_csv, "Points:0", "vel:1")
        ax.plot(x, v, color=color, linewidth=2, label=f"Re = {re}", zorder=3)
        ax.plot(
            ghia["y"], ghia[f"Re={re}"],
            linestyle="none", marker="o", markersize=6,
            markerfacecolor="none", markeredgecolor=color, markeredgewidth=1.3,
            zorder=4,
        )

    ax.set_xlabel("x")
    ax.set_ylabel("v")
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.6, 0.6)
    ax.set_title("v along horizontal centerline (y = 0.5)")
    style_axes(ax)
    ax.legend(loc="best", frameon=False)
    fig.text(0.5, 0.01, "line: present solver   o: Ghia et al. (1982)",
              ha="center", color=INK, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")
    if show:
        plt.show(block=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true",
                         help="save the PNGs but don't open plot windows")
    args = parser.parse_args()

    show = not args.no_show
    plot_u_mid(SCRIPT_DIR / "u_mid_comparison.png", show)
    plot_v_mid(SCRIPT_DIR / "v_mid_comparison.png", show)

    if show:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
