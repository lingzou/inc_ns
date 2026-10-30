#!/usr/bin/env python3
"""Report the backward-facing-step reattachment length, in meters, for each
Re case written by run_Re_sweep.py, and compare against Biswas et al.

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

Converged cases are additionally compared against the reference reattachment
lengths in Biswas_original.csv, which reports Lr/h normalized by the step
height (STEP_HEIGHT_M); this work's raw Lr is normalized the same way before
comparing. Written to comparison.csv and plotted as comparison.png.

Usage: python3 benchmark/backward_facing_step/report_reattachment.py
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file

THIS_DIR = Path(__file__).resolve().parent
# inc_exec writes "<input_base>_out..." alongside the input file (see
# cxx_main.C's output_prefix); run_Re_sweep.py writes the .i files into THIS_DIR.
OUTPUT_DIR = THIS_DIR
BISWAS_CSV = THIS_DIR / "Biswas_original.csv"
COMPARISON_CSV = THIS_DIR / "comparison.csv"
COMPARISON_PNG = THIS_DIR / "comparison.png"

RE_NUMBERS = [50, 100, 150, 200, 250, 300, 400]
TOL = 1e-6
CONVERGED_TOL = 1e-4  # final u_L2_error below this is "converged enough to trust Lr"
STEP_HEIGHT_M = 0.9423  # Biswas_original.csv reports Lr/h; normalize this work's Lr the same way

# dataviz reference palette (references/palette.md): categorical slots 1 & 2,
# validated for a 2-series light-mode scatter.
COLOR_THIS_WORK = "#2a78d6"  # slot 1, blue
COLOR_BISWAS = "#eb6834"  # slot 2, orange
COLOR_MUTED = "#898781"  # muted ink, for the parity line


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


def load_biswas(csv_path):
    """Return {Re: Lr/h} from Biswas_original.csv (plain 'Re, Lr/h' rows, no header)."""
    reference = {}
    with open(csv_path, newline="") as f:
        for row in csv.reader(f):
            if not row:
                continue
            reference[int(row[0])] = float(row[1])
    return reference


def write_comparison_csv(rows, csv_path):
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Re", "Lr/h_Biswas", "Lr/h_this_work", "error_abs", "error_pct [%]"])
        for Re_num, Lr_ref, Lr_work in rows:
            err_abs = Lr_work - Lr_ref
            err_pct = 100.0 * err_abs / Lr_ref
            writer.writerow([Re_num, f"{Lr_ref:.4f}", f"{Lr_work:.4f}", f"{err_abs:.4f}", f"{err_pct:.2f}"])


def plot_comparison(rows, png_path):
    Re_vals = [r[0] for r in rows]
    Lr_ref = [r[1] for r in rows]
    Lr_work = [r[2] for r in rows]

    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=300)
    ax.scatter(Re_vals, Lr_ref, s=64, marker="o", facecolors="none",
               edgecolors=COLOR_BISWAS, linewidths=1.5, label="Biswas et al.")
    ax.scatter(Re_vals, Lr_work, s=48, marker="x", color=COLOR_THIS_WORK,
               linewidths=1.5, label="This work")

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1)

    ax.set_xlabel("Re")
    ax.set_ylabel("Normalized reattachment length, $L_r/h$")
    ax.set_xlim(0, 500)
    ax.set_ylim(0, 10)
    ax.set_title("Backward-facing-step reattachment length")
    ax.grid(True, color="#e1e0d9", linewidth=0.8)
    ax.set_axisbelow(True)
    # for spine in ("top", "right"):
    #     ax.spines[spine].set_visible(False)
    ax.legend(frameon=True, loc="upper left")

    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def main():
    print(f"{'Re':>5}  {'Lr [m]':>10}  {'u_L2_error':>12}  status")
    converged_Lr = {}
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

        if converged and Lr is not None:
            converged_Lr[Re_num] = Lr

    if not BISWAS_CSV.is_file():
        print(f"\n{BISWAS_CSV.name} not found, skipping comparison")
        return

    biswas = load_biswas(BISWAS_CSV)
    rows = [(Re_num, biswas[Re_num], converged_Lr[Re_num] / STEP_HEIGHT_M)
            for Re_num in RE_NUMBERS if Re_num in biswas and Re_num in converged_Lr]

    if not rows:
        print(f"\nno converged cases overlap with {BISWAS_CSV.name}, skipping comparison")
        return

    write_comparison_csv(rows, COMPARISON_CSV)
    plot_comparison(rows, COMPARISON_PNG)

    print(f"\n{'Re':>5}  {'Biswas Lr/h':>11}  {'this work Lr/h':>14}  {'err [%]':>8}")
    for Re_num, Lr_ref, Lr_work in rows:
        err_pct = 100.0 * (Lr_work - Lr_ref) / Lr_ref
        print(f"{Re_num:>5}  {Lr_ref:>11.4f}  {Lr_work:>14.4f}  {err_pct:>8.2f}")
    print(f"\nwrote {COMPARISON_CSV.name} and {COMPARISON_PNG.name}")


if __name__ == "__main__":
    main()
