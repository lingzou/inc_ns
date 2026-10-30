#!/usr/bin/env python3
"""Extract lid-driven cavity centerline profiles for Ghia et al. (1982) validation.

Step 1: four converged Exodus results already sit next to this script, one per
Reynolds number: cavity_Re100_out.e, cavity_Re400_out.e, cavity_Re1000_out.e,
cavity_Re3200_out.e. With --newton, the four monolithic Newton counterparts
(cavity_Re<Re>_newton_out_newton.e, from run_Re_sweep.py --newton) are used
instead.

Step 2: for each one, shell out to ParaView's own Python (pvpython) to run
plot_lines.py, which samples u along the vertical centerline (x = 0.5) and v
along the horizontal centerline (y = 0.5) and writes them as CSV. plot_lines.py
always names its two outputs u_mid.csv/v_mid.csv, so each run happens in a
scratch directory and the results are moved here as u_mid_Re<Re>.csv /
v_mid_Re<Re>.csv (or u_mid_Re<Re>_newton.csv / v_mid_Re<Re>_newton.csv under
--newton). Skipped for a Re whose CSVs already exist, unless --force is given
-- pvpython startup dominates the runtime.

Usage:
    python3 benchmark/lid_driven/extract_lid_driven.py
    python3 benchmark/lid_driven/extract_lid_driven.py --newton
    python3 benchmark/lid_driven/extract_lid_driven.py --pvpython /path/to/pvpython
    python3 benchmark/lid_driven/extract_lid_driven.py --force
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_LINES_SCRIPT = SCRIPT_DIR / "plot_lines.py"

DEFAULT_PVPYTHON = "/Applications/ParaView-6.1.1.app/Contents/bin/pvpython"

REYNOLDS_NUMBERS = [100, 400, 1000, 3200]


def exodus_path(re, newton):
    if newton:
        return SCRIPT_DIR / f"cavity_Re{re}_newton_out_newton.e"
    return SCRIPT_DIR / f"cavity_Re{re}_out.e"


def sim_csv_paths(re, newton):
    suffix = "_newton" if newton else ""
    return SCRIPT_DIR / f"u_mid_Re{re}{suffix}.csv", SCRIPT_DIR / f"v_mid_Re{re}{suffix}.csv"


def extract_profiles(re, newton, pvpython, force):
    """Run plot_lines.py under pvpython for one Re, unless the CSVs already exist."""
    u_csv, v_csv = sim_csv_paths(re, newton)
    if u_csv.exists() and v_csv.exists() and not force:
        print(f"Re={re}: {u_csv.name}, {v_csv.name} already exist, skipping pvpython")
        return

    src = exodus_path(re, newton)
    if not src.exists():
        raise FileNotFoundError(f"missing Exodus output for Re={re}: {src}")

    with tempfile.TemporaryDirectory() as tmp_dir:
        cmd = [pvpython, str(PLOT_LINES_SCRIPT), str(src), tmp_dir]
        print(f"Re={re}: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        shutil.move(f"{tmp_dir}/u_mid.csv", u_csv)
        shutil.move(f"{tmp_dir}/v_mid.csv", v_csv)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pvpython", default=DEFAULT_PVPYTHON,
                         help=f"path to ParaView's pvpython (default: {DEFAULT_PVPYTHON})")
    parser.add_argument("--newton", action="store_true",
                         help="extract from the monolithic Newton Exodus outputs instead "
                              "of the segregated solver's")
    parser.add_argument("--force", action="store_true",
                         help="re-run pvpython even if the per-Re CSVs already exist")
    args = parser.parse_args()

    for re in REYNOLDS_NUMBERS:
        extract_profiles(re, args.newton, args.pvpython, args.force)


if __name__ == "__main__":
    sys.exit(main())
