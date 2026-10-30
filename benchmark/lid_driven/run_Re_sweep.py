#!/usr/bin/env python3
"""Generate lid-driven cavity input files for the Re sweep in cavity_base.i and run them.

Run from anywhere; it resolves paths relative to this file and executes
build/inc_exec from the repo root (required for the relative mesh path in
the input files to resolve).
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

# Re -> (visc, dt, n_timestep, output_interval), taken from the header comment
# in cavity_base.i.
CASES = {
    100:  dict(visc=0.01,      dt=0.01, n_timestep=3000,  output_interval=100),
    400:  dict(visc=0.0025,    dt=0.01, n_timestep=5000,  output_interval=200),
    1000: dict(visc=0.001,     dt=0.01, n_timestep=6000,  output_interval=200),
    3200: dict(visc=0.0003125, dt=0.01, n_timestep=10000, output_interval=200),
}

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
BASE_INPUT = THIS_DIR / "cavity_base.i"
INC_EXEC = REPO_ROOT / "build" / "inc_exec"

# Matrix-free Jacobian action + Schur-complement fieldsplit preconditioner,
# matching input/cavity_newton_mf_Schur.i. Appended verbatim to cavity_base.i
# for the --newton cases.
NEWTON_OPTIONS = """
############################################################
# Newton solver: matrix-free Jacobian action (-snes_mf_operator) with a
# Schur-complement fieldsplit preconditioner. See input/cavity_newton_mf_Schur.i
# for why each option below is set.
############################################################
-newton true
-snes_mf_operator

-ksp_type fgmres
-pc_type fieldsplit
-pc_fieldsplit_type schur
-pc_fieldsplit_schur_fact_type full
-pc_fieldsplit_0_fields 0,1
-pc_fieldsplit_1_fields 2

-fieldsplit_0_ksp_type preonly
-fieldsplit_0_pc_type ilu
-fieldsplit_0_pc_factor_levels 2
-fieldsplit_0_pc_factor_mat_ordering_type rcm

-fieldsplit_1_ksp_type preonly
-fieldsplit_1_pc_type gamg
-newton_schur_mp true
"""

# Per-Re overrides needed only under --newton, appended after NEWTON_OPTIONS.
# Re=3200: the solver's default initial guess (u=1, v=0 everywhere, see -u_init/-v_init
# in cxx_main.C) equals the lid velocity but violates no-slip at every other wall; at
# this case's low viscosity that mismatch makes the very first Newton linear solve fail
# (SNES diverged, reason -3, 0 iterations) before any nonlinear iteration completes.
# Starting from rest (u=0, v=0) avoids the issue -- confirmed to run cleanly. Re=100/400/
# 1000 converge fine with the default guess, so they don't need this override.
NEWTON_INIT_OVERRIDES = {
    3200: """
# Default initial guess (u=1, v=0) makes the first Newton linear solve fail at this
# Re's low viscosity -- start from rest instead.
-u_init 0.0
-v_init 0.0
""",
}


def sub_option(option, value, text):
    pattern = rf"^-{option} .*$"
    replacement = f"-{option} {value}"
    new_text, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if n != 1:
        raise RuntimeError(f"expected exactly one '-{option}' line in {BASE_INPUT}, found {n}")
    return new_text


def make_input(Re_num, params, newton, suffix):
    text = BASE_INPUT.read_text()
    for option in ("dt", "visc", "n_timestep", "output_interval"):
        text = sub_option(option, params[option], text)
    if newton:
        text = text.rstrip("\n") + "\n" + NEWTON_OPTIONS.lstrip("\n")
        if Re_num in NEWTON_INIT_OVERRIDES:
            text = text.rstrip("\n") + "\n" + NEWTON_INIT_OVERRIDES[Re_num].lstrip("\n")

    out_path = THIS_DIR / f"cavity_Re{Re_num}{suffix}.i"
    out_path.write_text(text)
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--newton", action="store_true",
        help="generate and run cavity_Re<Re>_newton.i cases using the matrix-free "
             "Schur Newton solver options (see NEWTON_OPTIONS) instead of the default solver",
    )
    args = parser.parse_args()
    suffix = "_newton" if args.newton else ""

    if not INC_EXEC.is_file():
        sys.exit(f"ERROR: {INC_EXEC} not found. Build first (see compile.sh).")

    input_files = {
        Re_num: make_input(Re_num, params, args.newton, suffix) for Re_num, params in CASES.items()
    }
    for Re_num, path in input_files.items():
        print(f"wrote {path.relative_to(REPO_ROOT)}")

    procs = {}
    log_files = {}
    for Re_num, path in input_files.items():
        log_path = THIS_DIR / f"cavity_Re{Re_num}{suffix}.log"
        log_file = open(log_path, "w")
        proc = subprocess.Popen(
            [str(INC_EXEC), "-i", str(path.relative_to(REPO_ROOT))],
            cwd=REPO_ROOT,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
        procs[Re_num] = proc
        log_files[Re_num] = log_file
        print(f"started Re={Re_num} (pid {proc.pid}), log: {log_path}")

    results = {}
    for Re_num, proc in procs.items():
        rc = proc.wait()
        log_files[Re_num].close()
        results[Re_num] = rc
        status = "OK" if rc == 0 else f"FAILED (rc={rc})"
        print(f"Re={Re_num}: {status}")

    ok = all(rc == 0 for rc in results.values())
    variant = " (Newton, matrix-free Schur)" if args.newton else ""
    summary = f"lid-driven cavity Re sweep{variant}: " + (
        "all 4 cases finished OK" if ok
        else "finished with failures: " + ", ".join(f"Re={Re_num}" for Re_num, rc in results.items() if rc != 0)
    )
    print(summary)

    if sys.platform == "darwin":
        subprocess.run(
            ["osascript", "-e", f'display notification "{summary}" with title "inc_ns benchmark"'],
            check=False,
        )


if __name__ == "__main__":
    main()
