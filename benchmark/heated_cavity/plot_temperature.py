#!/usr/bin/env python3
"""Temperature color map with contour lines for the Ra=1e3 heated cavity.

Hard-coded for cavity_buoyancy_Ra1e3_out.e: flat-shaded T color map (T is
stored elemental, piecewise-constant per triangle) overlaid with contour
lines from 0 to 1 in steps of 0.05.

Usage:
    python3 benchmark/heated_cavity/plot_temperature.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.tri import Triangulation
from scipy.io import netcdf_file

SCRIPT_DIR = Path(__file__).resolve().parent
##### Ra = 1e3
# EXODUS_FILE = SCRIPT_DIR / "cavity_buoyancy_Ra1e3_out.e"
# OUT_PNG = SCRIPT_DIR / "cavity_buoyancy_Ra1e3_temperature.png"

##### Ra = 1e4
# EXODUS_FILE = SCRIPT_DIR / "cavity_buoyancy_Ra1e4_out.e"
# OUT_PNG = SCRIPT_DIR / "cavity_buoyancy_Ra1e4_temperature.png"

##### Ra = 1e5
# EXODUS_FILE = SCRIPT_DIR / "cavity_buoyancy_Ra1e5_out.e"
# OUT_PNG = SCRIPT_DIR / "cavity_buoyancy_Ra1e5_temperature.png"

##### Ra = 1e6
EXODUS_FILE = SCRIPT_DIR / "cavity_buoyancy_Ra1e6_out.e"
OUT_PNG = SCRIPT_DIR / "cavity_buoyancy_Ra1e6_temperature.png"

LEVELS = np.arange(0.0, 1.0 + 1e-9, 0.05)

# Sequential blue ramp from the dataviz skill's reference palette.
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#2a78d6", "#0d366b"]
)
INK = "#0b0b0b"
MUTED = "#898781"


def name_str(char_array):
    return b"".join(char_array).decode("utf-8", errors="ignore").rstrip("\x00")


def elemental_to_nodal(values, triangles, num_nodes):
    """Average each triangle's constant value onto its three nodes.

    tricontour needs one scalar per node, but T is stored one-per-triangle.
    This average is only a guide for where the contour lines fall -- the
    tripcolor fill underneath still shows the real, piecewise-constant data.
    """
    sums = np.zeros(num_nodes)
    counts = np.zeros(num_nodes)
    flat_nodes = triangles.ravel()
    np.add.at(sums, flat_nodes, np.repeat(values, 3))
    np.add.at(counts, flat_nodes, 1)
    return sums / counts


def main():
    nc = netcdf_file(str(EXODUS_FILE), "r", mmap=False)
    x = nc.variables["coordx"][:]
    y = nc.variables["coordy"][:]
    conn = nc.variables["connect1"][:] - 1  # exodus connectivity is 1-indexed
    tri = Triangulation(x, y, conn)

    elem_names = [name_str(nc.variables["name_elem_var"][i])
                  for i in range(nc.variables["name_elem_var"].shape[0])]
    idx = elem_names.index("T") + 1
    T = nc.variables[f"vals_elem_var{idx}eb1"][-1]

    fig, ax = plt.subplots(figsize=(6, 6))

    mesh = ax.tripcolor(tri, facecolors=T, cmap=SEQUENTIAL,
                         vmin=0.0, vmax=1.0, shading="flat")

    T_nodal = elemental_to_nodal(T, tri.triangles, len(x))
    cs = ax.tricontour(tri, T_nodal, levels=LEVELS, colors=INK, linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=6, fmt="%.2f")

    ax.set_aspect("equal")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)

    cbar = fig.colorbar(mesh, ax=ax, shrink=0.85)
    cbar.set_label("T", color=INK)
    cbar.ax.tick_params(colors=MUTED)

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300)
    print(f"wrote {OUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
