#!/usr/bin/env python3
"""Report the average Nusselt number from a heated-cavity Ra sweep.

Reads the last row of each "cavity_buoyancy_Ra<Ra>_out_convergence.csv" (see
run_Ra_sweep.py, cxx_main.C's -nusselt_zones) written by build/inc_exec, prints
a table against de Vahl Davis's (1983) reference Nu_avg values, and plots it
on a log-x axis.

Nu_4 (hot wall) and Nu_2 (cold wall) should agree in magnitude at steady state
(energy balance); the reported "present solver" value is their average
magnitude, and the table's "hot/cold mismatch" column is a convergence sanity
check, not a modeling uncertainty.

Usage:
    python3 benchmark/heated_cavity/report_nusselt.py
    python3 benchmark/heated_cavity/report_nusselt.py --no-show
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent

RA_NUMBERS = [1e3, 1e4, 1e5, 1e6]

# de Vahl Davis (1983), "Natural convection of air in a square cavity: a bench
# mark numerical solution" -- the reference values every later paper on this
# case (including DEB08 Table 1) compares against.
DE_VAHL_DAVIS_NU = {1e3: 1.118, 1e4: 2.243, 1e5: 4.519, 1e6: 8.799}

# DEB08's own Table 1
DEB08_NU = {1e3: 1.160, 1e4: 2.247, 1e5: 4.613, 1e6: 9.106}

# Same visual language as benchmark/lid_driven/plot_lid_driven.py.
ACCENT = "#2a78d6"
REFERENCE = "#898781"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"
MUTED = "#898781"
INK = "#0b0b0b"


def ra_label(Ra):
    import math
    return f"1e{round(math.log10(Ra))}"


def convergence_csv(Ra):
    # inc_exec writes "<input_base>_out..." alongside the input file (see
    # cxx_main.C's output_prefix); run_Ra_sweep.py writes cavity_buoyancy_Ra<Ra>.i
    # next to this script.
    return SCRIPT_DIR / f"cavity_buoyancy_Ra{ra_label(Ra)}_out_convergence.csv"


def load_last_row(Ra):
    path = convergence_csv(Ra)
    if not path.exists():
        raise FileNotFoundError(f"missing {path} -- run run_Ra_sweep.py first")
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
    print(f"{'Ra':>8} {'Nu_hot':>9} {'Nu_cold':>9} {'mismatch':>9} {'present':>9} "
          f"{'deVahlDavis':>12} {'err%':>7} {'DEB08':>8} {'err%':>7} {'T_L2_error':>11}")
    for Ra in RA_NUMBERS:
        row = load_last_row(Ra)
        nu_hot = float(row["Nu_4"])
        nu_cold = -float(row["Nu_2"])  # cold wall is sign-flipped by convention
        present = 0.5 * (nu_hot + nu_cold)
        mismatch = abs(nu_hot - nu_cold) / present * 100.0
        ref = DE_VAHL_DAVIS_NU[Ra]
        err = (present - ref) / ref * 100.0
        deb08 = DEB08_NU.get(Ra)
        deb08_err = (present - deb08) / deb08 * 100.0 if deb08 is not None else float("nan")
        results[Ra] = present
        print(f"{Ra:>8.0e} {nu_hot:>9.4f} {nu_cold:>9.4f} {mismatch:>8.3f}% {present:>9.4f} "
              f"{ref:>12.3f} {err:>6.2f}% "
              f"{'' if deb08 is None else f'{deb08:.3f}':>8} "
              f"{'' if deb08 is None else f'{deb08_err:.2f}%':>7} {float(row['T_L2_error']):>11.2e}")

    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot(RA_NUMBERS, [results[Ra] for Ra in RA_NUMBERS], color=ACCENT, linewidth=2,
            marker="s", markersize=6, label="present solver", zorder=3)
    ax.plot(RA_NUMBERS, [DE_VAHL_DAVIS_NU[Ra] for Ra in RA_NUMBERS],
            linestyle="none", marker="o", markersize=9,
            markerfacecolor="none", markeredgecolor=REFERENCE, markeredgewidth=1.5,
            label="de Vahl Davis (1983)", zorder=4)
    deb08_ras = sorted(DEB08_NU)
    ax.plot(deb08_ras, [DEB08_NU[Ra] for Ra in deb08_ras],
            linestyle="none", marker="^", markersize=9,
            markerfacecolor="none", markeredgecolor=INK, markeredgewidth=1.5,
            label="DEB08 Table 1", zorder=5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Ra")
    ax.set_ylabel("Nu_avg")
    ax.set_title("Differentially heated cavity: average Nusselt number vs. Ra")
    ax.grid(True, which="both", color=GRIDLINE, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.legend(loc="best", frameon=False)

    out_png = SCRIPT_DIR / "nusselt_vs_Ra.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    print(f"\nwrote {out_png.relative_to(REPO_ROOT)}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
