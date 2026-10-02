# Lid-driven cavity benchmark

Validates the solver against Ghia et al. (1982) centerline profiles for
Re = 100, 400, 1000, 3200.

## Scripts

- `run_Re_sweep.py` — generates `cavity_Re<Re>.i` from `cavity_base.i` for
  each Re and runs `inc_exec` on all of them, producing `cavity_Re<Re>_out.e`.
- `extract_lid_driven.py` — runs pvpython over each `cavity_Re<Re>_out.e` to
  sample the u/v centerline profiles, writing `u_mid_Re<Re>.csv` /
  `v_mid_Re<Re>.csv`.
- `plot_lines.py` — pvpython helper invoked by `extract_lid_driven.py`; not
  run directly. Samples u/v along x = 0.5 / y = 0.5 for one Exodus file and
  writes `u_mid.csv` / `v_mid.csv`.
- `plot_lid_driven.py` — reads those CSVs plus `Ghia-u.csv`/`Ghia-v.csv` and
  produces `u_mid_comparison.png` / `v_mid_comparison.png`.
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
