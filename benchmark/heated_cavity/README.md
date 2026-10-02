# Differentially heated cavity benchmark

Validates the solver against de Vahl Davis (1983) average Nusselt numbers
for Ra = 1e3, 1e4, 1e5, 1e6 (Pr = 0.71).

## Scripts

- `run_Ra_sweep.py` — generates `cavity_buoyancy_Ra<Ra>.i` from
  `cavity_buoyancy_base.i` for each Ra (deriving `visc`/`diff_t` from Ra/Pr)
  and runs `build/inc_exec` on all of them in parallel, producing
  `cavity_buoyancy_Ra<Ra>_out.e` and `..._out_convergence.csv`.
- `report_nusselt.py` — reads the last row of each convergence CSV, prints a
  table comparing the solver's Nu_avg against de Vahl Davis (1983) and DEB08
  reference values, and writes `nusselt_vs_Ra.png`.
- `plot_temperature.py` — flat-shaded T color map with contour lines (0 to 1
  in steps of 0.05) for one `cavity_buoyancy_Ra<Ra>_out.e` file. Which Ra is
  plotted is hard-coded at the top of the script (edit `EXODUS_FILE` /
  `OUT_PNG` to switch).

## Usage

```
python3 run_Ra_sweep.py
python3 report_nusselt.py
python3 plot_temperature.py
```
