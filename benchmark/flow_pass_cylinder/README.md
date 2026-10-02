# Flow past a cylinder benchmark

Validates the solver's drag coefficient against Dennis & Chang (1970) for
the steady, symmetric regime: Re = 10, 20, 40.

## Scripts

- `run_Re_sweep.py` — generates `cylinder_Re<Re>.i` from `cylinder_base.i`
  for each Re (deriving `visc` from Re, with extra timesteps for Re = 20/40
  to reach convergence) and runs `build/inc_exec` on all of them in
  parallel, producing `cylinder_Re<Re>_out.e` and `..._out_convergence.csv`.
- `report_drag.py` — reads the last row of each convergence CSV, prints a
  table comparing the solver's C_D (split into pressure/friction parts)
  against Dennis & Chang (1970) and DEB08 reference values, and writes
  `drag_vs_Re.png`.
- `plot_field_contour.py` — contour-plots a field (default: pressure) from a
  `*_out.e` exodus file over the mesh triangulation, with labeled contour
  lines. Works for any benchmark case, not just this one.

## Usage

```
python3 run_Re_sweep.py
python3 report_drag.py
python3 plot_field_contour.py cylinder_Re10_out.e

python3 plot_field_contour.py cylinder_Re10_out.e --vmin -0.50 --vmax 0.70 --contours 23
python3 plot_field_contour.py cylinder_Re20_out.e --vmin -0.50 --vmax 0.60 --contours 21
python3 plot_field_contour.py cylinder_Re40_out.e --vmin -0.45 --vmax 0.55 --contours 19
```
