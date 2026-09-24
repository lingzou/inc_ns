#!/usr/bin/env python3
"""Report the backward-facing-step reattachment length, in meters, for each
Re case written by run_Re_sweep.py.

input/check_reattachment.py hardcodes the step location (X_STEP = 0.0) for
the *other* BFS meshes (see docs/to-do-list: item 9's mesh desync, items
2/3's "won't fix" on that script) -- wrong for this benchmark's mesh, where
the step is at x = 1. This instead detects the step location from each
case's own mesh: the bottom wall is the domain's lowest y, and the step
sits at the smallest x reached by nodes on that wall. No mesh-specific
constant needed, and it works unchanged if the mesh is ever swapped.

It also will not report a case whose run has not converged (the other gap
check_reattachment.py's own docstring calls out): CONVERGED_TOL below is a
practical threshold on the final per-step u_L2_error, not a claim that the
whole field is converged to machine precision.

Usage: python3 benchmark/backward_facing_step/report_reattachment.py
"""
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

THIS_DIR = Path(__file__).resolve().parent
# inc_exec writes "<input_base>_out..." alongside the input file (see
# cxx_main.C's output_prefix); run_Re_sweep.py writes the .i files into THIS_DIR.
OUTPUT_DIR = THIS_DIR

RE_NUMBERS = [50, 100, 150, 200, 250, 300, 400]
TOL = 1e-6
CONVERGED_TOL = 1e-4  # final u_L2_error below this is "converged enough to trust Lr"


def reattachment_length_m(e_path):
    f = netcdf_file(e_path, "r", mmap=False)
    x = f.variables["coordx"][:]
    y = f.variables["coordy"][:]
    conn = f.variables["connect1"][:] - 1  # 1-based in Exodus

    y_wall = y.min()
    on_wall_node = np.abs(y - y_wall) < TOL
    x_step = x[on_wall_node].min()

    on_wall = np.abs(y[conn] - y_wall) < TOL
    wall_cells = on_wall.sum(axis=1) == 2  # cells with an edge on the wall
    cx = x[conn].mean(axis=1)
    downstream = wall_cells & (cx > x_step)

    names = [b"".join(r).decode().strip() for r in f.variables["name_elem_var"][:]]
    u_idx = names.index("vel_x") + 1
    last = f.variables["time_whole"].shape[0]
    u_all = f.variables[f"vals_elem_var{u_idx}eb1"][last - 1]

    order = np.argsort(cx[downstream])
    xs = cx[downstream][order]
    us = u_all[downstream][order]

    for i in range(len(us) - 1):
        if us[i] < 0 and us[i + 1] > 0:
            frac = -us[i] / (us[i + 1] - us[i])
            x_r = xs[i] + frac * (xs[i + 1] - xs[i])
            return x_r - x_step

    return None


def final_u_l2_error(csv_path):
    if not csv_path.is_file():
        return None
    last_line = csv_path.read_text().strip().splitlines()[-1]
    return float(last_line.split(",")[1])


def main():
    print(f"{'Re':>5}  {'Lr [m]':>10}  {'u_L2_error':>12}  status")
    for Re_num in RE_NUMBERS:
        e_path = OUTPUT_DIR / f"backward_facing_step_Re{Re_num}_out.e"
        csv_path = OUTPUT_DIR / f"backward_facing_step_Re{Re_num}_out_convergence.csv"

        if not e_path.is_file():
            print(f"{Re_num:>5}  {'--':>10}  {'--':>12}  missing output ({e_path.name})")
            continue

        Lr = reattachment_length_m(e_path)
        err = final_u_l2_error(csv_path)
        converged = err is not None and err < CONVERGED_TOL
        Lr_str = f"{Lr:.4f}" if Lr is not None else "n/a"
        err_str = f"{err:.3g}" if err is not None else "?"
        status = "ok" if converged else "NOT CONVERGED"
        print(f"{Re_num:>5}  {Lr_str:>10}  {err_str:>12}  {status}")


if __name__ == "__main__":
    main()
