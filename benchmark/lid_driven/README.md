# Lid-driven cavity benchmark

Validates the solver against Ghia et al. (1982) centerline profiles for
Re = 100, 400, 1000, 3200.

## Scripts

- `run_Re_sweep.py` — generates `cavity_Re<Re>.i` from `cavity_base.i` for
  each Re and runs `inc_exec` on all of them, producing `cavity_Re<Re>_out.e`.
  With `--newton`, generates `cavity_Re<Re>_newton.i` instead, adding the
  matrix-free Schur-complement Newton solver options (see
  `input/cavity_newton_mf_Schur.i`). The Re=3200 case also gets `-u_init 0.0
  -v_init 0.0`: the solver's default initial guess (u=1, v=0 everywhere)
  violates no-slip walls badly enough at this Re's low viscosity to fail the
  first Newton linear solve, so that case starts from rest instead.
- `extract_lid_driven.py` — runs pvpython over each `cavity_Re<Re>_out.e` to
  sample the u/v centerline profiles, writing `u_mid_Re<Re>.csv` /
  `v_mid_Re<Re>.csv`. With `--newton`, reads the Newton Exodus outputs instead
  (`cavity_Re<Re>_newton_out_newton.e`) and writes `u_mid_Re<Re>_newton.csv` /
  `v_mid_Re<Re>_newton.csv`.
- `plot_lines.py` — pvpython helper invoked by `extract_lid_driven.py`; not
  run directly. Samples u/v along x = 0.5 / y = 0.5 for one Exodus file and
  writes `u_mid.csv` / `v_mid.csv`.
- `plot_lid_driven.py` — reads those CSVs plus `Ghia-u.csv`/`Ghia-v.csv` and
  produces `u_mid_comparison.png` / `v_mid_comparison.png`. Also reads the
  `_newton` CSVs (run `extract_lid_driven.py --newton` first) and produces
  `u_mid_solver_comparison.png` / `v_mid_solver_comparison.png`, overlaying the
  segregated (SIMPLE) and monolithic (Newton) profiles directly -- solid vs.
  dashed line for the same Re color -- to check the two solvers agree.
- `plot_streamfunction.py` — labeled StreamFunction contour plot tracing the
  cavity's vortices for one `cavity_Re<Re>_out.e` file. Which Re and contour
  levels are plotted is hard-coded at the top of the script (edit
  `EXODUS_FILE` / `OUT_PNG` / `LEVELS` to switch).

## Usage

```
python3 run_Re_sweep.py
python3 extract_lid_driven.py
python3 plot_lid_driven.py
python3 plot_streamfunction.py
```

For the Newton matrix-free Schur solver instead, and to compare it against the
segregated solver above:

```
python3 run_Re_sweep.py --newton
python3 extract_lid_driven.py --newton
python3 plot_lid_driven.py
```
