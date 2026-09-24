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

## Usage

```
python3 run_Re_sweep.py
python3 report_drag.py
```

(Other scripts in this folder — `run_blockage_study.py`,
`run_refinement_study.py`, `run_upstream_study.py`, and their
`report_*.py` counterparts — are separate mesh-sensitivity studies, not
part of this basic sweep.)
