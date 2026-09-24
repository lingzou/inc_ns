#!/usr/bin/env python3
"""Report the cylinder drag coefficient from a flow-past-a-cylinder Re sweep.

Reads the last row of each "cylinder_Re<Re>_out_convergence.csv" (see
run_Re_sweep.py, cxx_main.C's -force_zones) written by build/inc_exec, prints a
table against the standard steady-regime reference values, and plots C_D vs Re
on log-log axes.

C_D is split into its pressure (Cp) and friction (Cf) parts, which the solver
reports separately; the pressure share rising with Re is the physical signature
of the case and a cheap sanity check on the split.

Two caveats the table's "err%" columns do not capture:

  * C_D here runs high -- ~+13% against DEB08 at Re = 40 -- a known, tracked
    discrepancy (docs/to-do-list-09-22-2026.md item 13), with mesh density the
    working hypothesis. Read a change in err% across a solver change as signal;
    read its absolute size as that open item, not as a new bug.
  * The reference values are for an *unbounded* cylinder, while this mesh is a
    13D-tall channel (7.7% blockage), which raises C_D on its own. Part of the
    offset above is physics, not discretization error.

C_L must be zero at these Re by symmetry, so the printed value sizes the
first-order surface-force integration error rather than any lift.

Usage:
    python3 benchmark/flow_pass_cylinder/report_drag.py
    python3 benchmark/flow_pass_cylinder/report_drag.py --no-show
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent

RE_NUMBERS = [10, 20, 40]

# Dennis & Chang (1970), "Numerical solutions for steady flow past a circular
# cylinder at Reynolds numbers up to 100" -- the steady-regime reference set
# cited by DEB08 and by most later papers on this case. Unbounded flow.
DENNIS_CHANG_CD = {10: 2.846, 20: 2.045, 40: 1.522}
DENNIS_CHANG_Cf = {10: 1.246, 20: 0.812, 40: 0.522}
DENNIS_CHANG_Cp = {10: 1.600, 20: 1.233, 40: 0.998}

# DEB08's own Table 2
DEB08_CD       = {10: 2.968, 20: 2.146, 40: 1.599}
DEB08_CHANG_Cf = {10: 1.288, 20: 0.838, 40: 0.537}
DEB08_CHANG_Cp = {10: 1.680, 20: 1.308, 40: 1.062}

# Same visual language as benchmark/lid_driven/plot_lid_driven.py.
ACCENT = "#2a78d6"
REFERENCE = "#898781"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"
MUTED = "#898781"
INK = "#0b0b0b"


def convergence_csv(Re_num):
    # inc_exec writes "<input_base>_out..." alongside the input file (see
    # cxx_main.C's output_prefix); run_Re_sweep.py writes cylinder_Re<Re>.i next
    # to this script.
    return SCRIPT_DIR / f"cylinder_Re{Re_num}_out_convergence.csv"


def load_last_row(Re_num):
    path = convergence_csv(Re_num)
    if not path.exists():
        raise FileNotFoundError(f"missing {path} -- run run_Re_sweep.py first")
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"{path} has no data rows")
    return rows[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true",
                        help="save the PNG but don't open a plot window")
    args = parser.parse_args()

    results = {}
    print(f"{'Re':>4} {'Cp':>8} {'Cf':>8} {'Cp/Cd':>7} {'present':>8} "
          f"{'DennisChang':>12} {'err%':>7} {'DEB08':>8} {'err%':>7} "
          f"{'Cl':>10} {'u_L2_error':>11}")
    for Re_num in RE_NUMBERS:
        row = load_last_row(Re_num)
        cp = float(row["Cp"])
        cf = float(row["Cf"])
        present = float(row["Cd"])
        ref = DENNIS_CHANG_CD[Re_num]
        err = (present - ref) / ref * 100.0
        deb08 = DEB08_CD.get(Re_num)
        deb08_err = (present - deb08) / deb08 * 100.0 if deb08 is not None else None
        results[Re_num] = present
        print(f"{Re_num:>4} {cp:>8.4f} {cf:>8.4f} {cp / present:>7.3f} {present:>8.4f} "
              f"{ref:>12.3f} {err:>6.2f}% "
              f"{'' if deb08 is None else f'{deb08:.3f}':>8} "
              f"{'' if deb08_err is None else f'{deb08_err:.2f}%':>7} "
              f"{float(row['Cl']):>+10.2e} {float(row['u_L2_error']):>11.2e}")

    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot(RE_NUMBERS, [results[Re_num] for Re_num in RE_NUMBERS], color=ACCENT,
            linewidth=2, marker="s", markersize=6, label="present solver", zorder=3)
    ax.plot(RE_NUMBERS, [DENNIS_CHANG_CD[Re_num] for Re_num in RE_NUMBERS],
            linestyle="none", marker="o", markersize=9,
            markerfacecolor="none", markeredgecolor=REFERENCE, markeredgewidth=1.5,
            label="Dennis & Chang (1970)", zorder=4)
    deb08_res = sorted(DEB08_CD)
    ax.plot(deb08_res, [DEB08_CD[Re_num] for Re_num in deb08_res],
            linestyle="none", marker="^", markersize=9,
            markerfacecolor="none", markeredgecolor=INK, markeredgewidth=1.5,
            label="DEB08 Table 2", zorder=5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    # Only three Re, so label them explicitly and silence the log locator's
    # minor labels, which would otherwise print a stray "3 x 10^1" between them.
    ax.set_xticks(RE_NUMBERS)
    ax.set_xticklabels([str(Re_num) for Re_num in RE_NUMBERS])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Re")
    ax.set_ylabel("C_D")
    ax.set_title("Steady flow past a cylinder: drag coefficient vs. Re")
    ax.grid(True, which="both", color=GRIDLINE, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.legend(loc="best", frameon=False)

    out_png = SCRIPT_DIR / "drag_vs_Re.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    print(f"\nwrote {out_png.relative_to(REPO_ROOT)}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
