#!/usr/bin/env python3
"""Generate lid-driven cavity input files for the Re sweep in cavity_base.i and run them.

Run from anywhere; it resolves paths relative to this file and executes
build/inc_exec from the repo root (required for the relative mesh path in
the input files to resolve).
"""
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


def sub_option(option, value, text):
    pattern = rf"^-{option} .*$"
    replacement = f"-{option} {value}"
    new_text, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if n != 1:
        raise RuntimeError(f"expected exactly one '-{option}' line in {BASE_INPUT}, found {n}")
    return new_text


def make_input(Re_num, params):
    text = BASE_INPUT.read_text()
    for option in ("dt", "visc", "n_timestep", "output_interval"):
        text = sub_option(option, params[option], text)

    out_path = THIS_DIR / f"cavity_Re{Re_num}.i"
    out_path.write_text(text)
    return out_path


def main():
    if not INC_EXEC.is_file():
        sys.exit(f"ERROR: {INC_EXEC} not found. Build first (see compile.sh).")

    input_files = {Re_num: make_input(Re_num, params) for Re_num, params in CASES.items()}
    for Re_num, path in input_files.items():
        print(f"wrote {path.relative_to(REPO_ROOT)}")

    procs = {}
    log_files = {}
    for Re_num, path in input_files.items():
        log_path = THIS_DIR / f"cavity_Re{Re_num}.log"
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
    summary = "lid-driven cavity Re sweep: " + (
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
