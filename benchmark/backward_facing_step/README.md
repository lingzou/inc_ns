# Backward-facing-step benchmark

Validates the solver's reattachment length against Armaly et al. for
Re = 50, 100, 150, 200, 250, 300, 400 (fixed inlet velocity, Re varied via
viscosity).

## Scripts

- `run_Re_sweep.py` — generates `backward_facing_step_Re<Re>.i` from
  `backward_facing_step_base.i` for each Re (varying `visc`, with extra
  timesteps for Re = 250/300/400 to reach convergence) and runs
  `build/inc_exec` on all of them in parallel, producing
  `backward_facing_step_Re<Re>_out.e` and `..._out_convergence.csv`.
- `report_reattachment.py` — reads each `_out.e` file, detects the step
  location from the mesh, and computes the reattachment length Lr from the
  wall-adjacent u velocity sign change. Also checks the final u_L2_error
  against `CONVERGED_TOL` and flags any case that hasn't converged. For
  converged cases, compares Lr against the reference values in `Biswas.csv`,
  writing `comparison.csv` (Biswas, this work, absolute/percent error) and a
  scatter plot `comparison.png`.

## Usage

```
python3 run_Re_sweep.py
python3 report_reattachment.py
```

`report_reattachment.py` requires `scipy`
