#!/usr/bin/env python3
"""Plot lid-driven cavity solver profiles against Ghia et al. (1982).

Reads the per-Re centerline CSVs (u_mid_Re<Re>.csv / v_mid_Re<Re>.csv)
produced by extract_lid_driven.py, overlays them on the Ghia-u.csv / Ghia-v.csv
reference tables, and produces two comparison plots (u along y, v along x),
each shown on screen and saved as a 300 dpi PNG.

It also reads the monolithic-Newton counterparts (u_mid_Re<Re>_newton.csv /
v_mid_Re<Re>_newton.csv, from `extract_lid_driven.py --newton`) and produces two
more plots overlaying the segregated (SIMPLE) and monolithic (Newton) profiles
directly, one line style per solver, to check the two solvers agree.

Finally, it reads u_mid_Re<Re>_t0.2[_newton].csv / v_mid_Re<Re>_t0.2[_newton].csv,
a one-off early-transient snapshot (20 steps at dt=0.01 -> t=0.2 s, both solvers
started from rest via -u_init 0 -v_init 0 for an apples-to-apples IC -- see
cavity_Re<Re>_t0.2[_newton].i) and produces two more segregated-vs-monolithic
plots at that fixed early time, rather than at steady state. These CSVs are not
produced by extract_lid_driven.py; generate them the same way (plot_lines.py
under pvpython) if regenerating from scratch.

Usage:
    python3 benchmark/lid_driven/extract_lid_driven.py           # first, to generate the CSVs
    python3 benchmark/lid_driven/extract_lid_driven.py --newton  # and the Newton CSVs
    python3 benchmark/lid_driven/plot_lid_driven.py
    python3 benchmark/lid_driven/plot_lid_driven.py --no-show
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

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


def sim_csv_paths(re, newton=False, t02=False):
    suffix = ("_t0.2" if t02 else "") + ("_newton" if newton else "")
    return SCRIPT_DIR / f"u_mid_Re{re}{suffix}.csv", SCRIPT_DIR / f"v_mid_Re{re}{suffix}.csv"


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


# Solver variant encoding for the segregated-vs-monolithic comparison plots. Re
# still owns color on the segregated line (identity, fixed order, matching every
# other plot in this file). The Newton line breaks that rule on purpose: it's
# drawn in one fixed, non-Re color and on top (higher zorder) so that where the
# two solutions overlap exactly -- the expected result -- the dashed NEWTON_COLOR
# line is visibly sitting on top of (and hiding) the solid Re-colored one, rather
# than the two being indistinguishable.
SOLVER_LINESTYLES = {
    "segregated": "-",
    "newton": "--",
}
NEWTON_COLOR = "#000000"
SOLVER_CAPTION = "solid (color): segregated (SIMPLE)   dashed (black, drawn on top): monolithic (Newton)"


def solver_comparison_legend_handles():
    handles = [
        Line2D([0], [0], color=COLORS[re], linewidth=2, linestyle=SOLVER_LINESTYLES["segregated"],
               label=f"Re = {re}")
        for re in REYNOLDS_NUMBERS
    ]
    handles.append(
        Line2D([0], [0], color=NEWTON_COLOR, linewidth=2, linestyle=SOLVER_LINESTYLES["newton"],
               label="Monolithic (Newton)")
    )
    return handles


def plot_u_mid_solver_comparison(out_png, show):
    fig, ax = plt.subplots(figsize=(6, 6))
    for re in REYNOLDS_NUMBERS:
        color = COLORS[re]
        u_csv, _ = sim_csv_paths(re, newton=False)
        y, u = load_sim_profile(u_csv, "Points:1", "vel:0")
        ax.plot(u, y, color=color, linewidth=2, linestyle=SOLVER_LINESTYLES["segregated"], zorder=3)

        u_csv_newton, _ = sim_csv_paths(re, newton=True)
        y_n, u_n = load_sim_profile(u_csv_newton, "Points:1", "vel:0")
        ax.plot(u_n, y_n, color=NEWTON_COLOR, linewidth=2, linestyle=SOLVER_LINESTYLES["newton"], zorder=4)

    ax.set_xlabel("u")
    ax.set_ylabel("y")
    ax.set_xlim(-0.6, 1.2)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title("u along vertical centerline (x = 0.5): segregated vs monolithic")
    style_axes(ax)
    ax.legend(handles=solver_comparison_legend_handles(), loc="best", frameon=False)
    fig.text(0.5, 0.01, SOLVER_CAPTION, ha="center", color=INK, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")
    if show:
        plt.show(block=False)


def plot_v_mid_solver_comparison(out_png, show):
    fig, ax = plt.subplots(figsize=(6, 6))
    for re in REYNOLDS_NUMBERS:
        color = COLORS[re]
        _, v_csv = sim_csv_paths(re, newton=False)
        x, v = load_sim_profile(v_csv, "Points:0", "vel:1")
        ax.plot(x, v, color=color, linewidth=2, linestyle=SOLVER_LINESTYLES["segregated"], zorder=3)

        _, v_csv_newton = sim_csv_paths(re, newton=True)
        x_n, v_n = load_sim_profile(v_csv_newton, "Points:0", "vel:1")
        ax.plot(x_n, v_n, color=NEWTON_COLOR, linewidth=2, linestyle=SOLVER_LINESTYLES["newton"], zorder=4)

    ax.set_xlabel("x")
    ax.set_ylabel("v")
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.6, 0.6)
    ax.set_title("v along horizontal centerline (y = 0.5): segregated vs monolithic")
    style_axes(ax)
    ax.legend(handles=solver_comparison_legend_handles(), loc="best", frameon=False)
    fig.text(0.5, 0.01, SOLVER_CAPTION, ha="center", color=INK, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")
    if show:
        plt.show(block=False)


def plot_u_mid_t02_comparison(out_png, show):
    fig, ax = plt.subplots(figsize=(6, 6))
    for re in REYNOLDS_NUMBERS:
        color = COLORS[re]
        u_csv, _ = sim_csv_paths(re, newton=False, t02=True)
        y, u = load_sim_profile(u_csv, "Points:1", "vel:0")
        ax.plot(u, y, color=color, linewidth=2, linestyle=SOLVER_LINESTYLES["segregated"], zorder=3)

        u_csv_newton, _ = sim_csv_paths(re, newton=True, t02=True)
        y_n, u_n = load_sim_profile(u_csv_newton, "Points:1", "vel:0")
        ax.plot(u_n, y_n, color=NEWTON_COLOR, linewidth=2, linestyle=SOLVER_LINESTYLES["newton"], zorder=4)

    ax.set_xlabel("u")
    ax.set_ylabel("y")
    ax.set_title("u along vertical centerline (x = 0.5), t = 0.2 s")
    style_axes(ax)
    ax.legend(handles=solver_comparison_legend_handles(), loc="best", frameon=False)
    fig.text(0.5, 0.01, SOLVER_CAPTION, ha="center", color=INK, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")
    if show:
        plt.show(block=False)


def plot_v_mid_t02_comparison(out_png, show):
    fig, ax = plt.subplots(figsize=(6, 6))
    for re in REYNOLDS_NUMBERS:
        color = COLORS[re]
        _, v_csv = sim_csv_paths(re, newton=False, t02=True)
        x, v = load_sim_profile(v_csv, "Points:0", "vel:1")
        ax.plot(x, v, color=color, linewidth=2, linestyle=SOLVER_LINESTYLES["segregated"], zorder=3)

        _, v_csv_newton = sim_csv_paths(re, newton=True, t02=True)
        x_n, v_n = load_sim_profile(v_csv_newton, "Points:0", "vel:1")
        ax.plot(x_n, v_n, color=NEWTON_COLOR, linewidth=2, linestyle=SOLVER_LINESTYLES["newton"], zorder=4)

    ax.set_xlabel("x")
    ax.set_ylabel("v")
    ax.set_title("v along horizontal centerline (y = 0.5), t = 0.2 s")
    style_axes(ax)
    ax.legend(handles=solver_comparison_legend_handles(), loc="best", frameon=False)
    fig.text(0.5, 0.01, SOLVER_CAPTION, ha="center", color=INK, fontsize=9)
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
    plot_u_mid_solver_comparison(SCRIPT_DIR / "u_mid_solver_comparison.png", show)
    plot_v_mid_solver_comparison(SCRIPT_DIR / "v_mid_solver_comparison.png", show)
    plot_u_mid_t02_comparison(SCRIPT_DIR / "u_mid_t0.2_solver_comparison.png", show)
    plot_v_mid_t02_comparison(SCRIPT_DIR / "v_mid_t0.2_solver_comparison.png", show)

    if show:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
