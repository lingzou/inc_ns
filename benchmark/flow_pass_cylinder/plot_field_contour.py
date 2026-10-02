#!/usr/bin/env python3
"""Contour-plot a field from an inc_ns exodus output file.

Reads an "*_out.e" file with scipy.io.netcdf_file (exodus is classic NetCDF,
so no exodus/netCDF4 package is needed) and renders one field over the mesh
triangulation. Works for any of the benchmark cases (flow_pass_cylinder,
backward_facing_step, heated_cavity, ...) since they all share the same
variable layout: vel_x, vel_y, pressure, T as elemental fields, plus
StreamFunction as a nodal field.

Usage:
    python3 benchmark/plot_field_contour.py benchmark/flow_pass_cylinder/cylinder_Re10_out.e
    python3 benchmark/plot_field_contour.py path/to/out.e --field vel_x --timestep 5
    python3 benchmark/plot_field_contour.py path/to/out.e --no-show
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.tri import Triangulation
from scipy.io import netcdf_file

# Diverging pair (blue <-> red, gray midpoint) and sequential blue ramp from
# the dataviz skill's reference palette.
DIVERGING = LinearSegmentedColormap.from_list(
    "blue_gray_red", ["#2a78d6", "#f0efec", "#e34948"]
)
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#2a78d6", "#0d366b"]
)
INK = "#0b0b0b"
MUTED = "#898781"


def name_str(char_array):
    return b"".join(char_array).decode("utf-8", errors="ignore").rstrip("\x00")


def elemental_to_nodal(values, triangles, num_nodes):
    """Average each triangle's constant value onto its three nodes.

    tricontour needs one scalar per node, but pressure/vel_x/... are stored
    one-per-triangle. This average is only a guide for where the contour
    lines fall -- the tripcolor fill underneath still shows the real,
    piecewise-constant data.
    """
    sums = np.zeros(num_nodes)
    counts = np.zeros(num_nodes)
    flat_nodes = triangles.ravel()
    np.add.at(sums, flat_nodes, np.repeat(values, 3))
    np.add.at(counts, flat_nodes, 1)
    return sums / counts


def load_field(nc, field, timestep):
    """Return (values, is_nodal) for `field` at `timestep` (-1 = last)."""
    elem_names = [name_str(nc.variables["name_elem_var"][i])
                  for i in range(nc.variables["name_elem_var"].shape[0])]
    nod_names = [name_str(nc.variables["name_nod_var"][i])
                 for i in range(nc.variables["name_nod_var"].shape[0])]

    if field in elem_names:
        idx = elem_names.index(field) + 1
        return nc.variables[f"vals_elem_var{idx}eb1"][timestep], False
    if field in nod_names:
        idx = nod_names.index(field) + 1
        return nc.variables[f"vals_nod_var{idx}"][timestep], True

    raise ValueError(
        f"field '{field}' not found; available: {elem_names + nod_names}"
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("exodus_file", type=Path, help="path to an *_out.e file")
    parser.add_argument("--field", default="pressure",
                         help="field name (default: pressure)")
    parser.add_argument("--timestep", type=int, default=-1,
                         help="time step index, -1 = last (default: -1)")
    parser.add_argument("--out", type=Path, default=None,
                         help="output PNG path (default: alongside the input file)")
    parser.add_argument("--contours", type=int, default=10,
                         help="number of labeled contour lines, 0 to disable (default: 10)")
    parser.add_argument("--vmin", type=float, default=None,
                         help="color/contour scale minimum (default: data minimum)")
    parser.add_argument("--vmax", type=float, default=None,
                         help="color/contour scale maximum (default: data maximum)")
    parser.add_argument("--no-show", action="store_true",
                         help="save the PNG but don't open a plot window")
    args = parser.parse_args()

    if not args.exodus_file.exists():
        parser.error(f"no such file: {args.exodus_file}")

    nc = netcdf_file(str(args.exodus_file), "r", mmap=False)
    x = nc.variables["coordx"][:]
    y = nc.variables["coordy"][:]
    conn = nc.variables["connect1"][:] - 1  # exodus connectivity is 1-indexed
    tri = Triangulation(x, y, conn)

    times = nc.variables["time_whole"][:]
    step = args.timestep if args.timestep >= 0 else len(times) + args.timestep
    values, is_nodal = load_field(nc, args.field, step)
    time = times[step]

    fig, ax = plt.subplots(figsize=(8, 6))

    vmin = args.vmin if args.vmin is not None else float(values.min())
    vmax = args.vmax if args.vmax is not None else float(values.max())
    if vmin < 0.0 < vmax:
        cmap, norm = DIVERGING, TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
    else:
        cmap, norm = SEQUENTIAL, plt.Normalize(vmin=vmin, vmax=vmax)

    if is_nodal:
        mesh = ax.tricontourf(tri, values, levels=40, cmap=cmap, norm=norm)
    else:
        # Elemental data is piecewise-constant per triangle, so a flat-shaded
        # tripcolor plots what was actually solved for; tricontourf would
        # need a node average and imply smoothness that isn't there.
        mesh = ax.tripcolor(tri, facecolors=values, cmap=cmap, norm=norm,
                             shading="flat")

    if args.contours > 0:
        nodal_values = values if is_nodal else elemental_to_nodal(
            values, tri.triangles, len(x)
        )

        levels = np.linspace(vmin, vmax, args.contours + 2)
        # avoid things like '-5.55111512e-17' which should just be 0
        levels = np.round(levels, decimals=10) + 0.0

        print(levels)

        cs = ax.tricontour(tri, nodal_values, levels=levels,
                            colors="black", linewidths=0.6, alpha=0.6)
        ax.clabel(cs, inline=False, fontsize=6, fmt="%.3g")

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(20.0, 30.0)
    ax.set_ylim(7.0, 18.0)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    # I don't need a title
    # ax.set_title(f"{args.field}  ({args.exodus_file.name}, t = {time:g})", color=INK)

    cbar = fig.colorbar(mesh, ax=ax, shrink=0.85)
    cbar.set_label(args.field, color=INK)
    cbar.ax.tick_params(colors=MUTED)

    fig.tight_layout()

    out_png = args.out or args.exodus_file.with_name(
        f"{args.exodus_file.stem}_{args.field}.png"
    )
    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
