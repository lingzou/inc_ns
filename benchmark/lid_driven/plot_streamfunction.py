#!/usr/bin/env python3
"""Contour plot of the primary-vortex streamlines in the Re=100 lid-driven cavity.

Hard-coded for cavity_Re100_out.e: draws labeled StreamFunction contour
lines at -0.01, -0.03, -0.05, -0.1, the classic levels used to trace the
primary vortex in Ghia et al. (1982)-style cavity plots.

Usage:
    python3 benchmark/lid_driven/plot_streamfunction.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.io import netcdf_file

SCRIPT_DIR = Path(__file__).resolve().parent

###### Re = 100
# EXODUS_FILE = SCRIPT_DIR / "cavity_Re100_out.e"
# OUT_PNG = SCRIPT_DIR / "cavity_Re100_streamfunction.png"
# LEVELS = [-0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -1e-4, -1e-5, 1e-6, 1e-5]

###### Re = 400
# EXODUS_FILE = SCRIPT_DIR / "cavity_Re400_out.e"
# OUT_PNG = SCRIPT_DIR / "cavity_Re400_streamfunction.png"
# LEVELS = [-0.11, -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -1e-4, -1e-5, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3]

###### Re = 1000
# EXODUS_FILE = SCRIPT_DIR / "cavity_Re1000_out.e"
# OUT_PNG = SCRIPT_DIR / "cavity_Re1000_streamfunction.png"
# LEVELS = [-0.1175, -0.115, -0.11, -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -1e-4, -1e-5, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 1.5e-3]

###### Re = 3200
###### The results is not well resolved as in Ghia's paper, but looks ok.
EXODUS_FILE = SCRIPT_DIR / "cavity_Re3200_out.e"
OUT_PNG = SCRIPT_DIR / "cavity_Re3200_streamfunction.png"
LEVELS = [-0.108, -0.105, -0.1, -0.09, -0.07, -0.05, -0.03, -0.01, -1e-4, 1e-5, 1e-4, 5e-4, 1e-3, 2e-3, 3e-3]

INK = "#0b0b0b"
MUTED = "#898781"


def name_str(char_array):
    return b"".join(char_array).decode("utf-8", errors="ignore").rstrip("\x00")


def main():
    nc = netcdf_file(str(EXODUS_FILE), "r", mmap=False)
    x = nc.variables["coordx"][:]
    y = nc.variables["coordy"][:]
    conn = nc.variables["connect1"][:] - 1  # exodus connectivity is 1-indexed
    tri = Triangulation(x, y, conn)

    nod_names = [name_str(nc.variables["name_nod_var"][i])
                 for i in range(nc.variables["name_nod_var"].shape[0])]
    idx = nod_names.index("StreamFunction") + 1
    psi = nc.variables[f"vals_nod_var{idx}"][-1]

    fig, ax = plt.subplots(figsize=(6, 6))

    cs = ax.tricontour(tri, psi, levels=LEVELS, colors=INK, linewidths=1.0)
    ax.clabel(cs, inline=True, fontsize=6, fmt="%.2e")

    ax.set_aspect("equal")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    # ax.set_title("Lid-driven cavity, Re = 100: StreamFunction contours")
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300)
    print(f"wrote {OUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
