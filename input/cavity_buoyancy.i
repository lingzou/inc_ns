-problem cavity_buoyancy
-mesh input/cavity.msh
-dt 0.05
-n_timestep 300 # 3000 for st-st
-output_interval 100

# Ra = 1e3, Pr = 0.71 differentially heated cavity (DEB08 section 4.5 / de Vahl Davis
# benchmark), nondimensionalized with -visc = sqrt(Pr/Ra), -diff_t = 1/sqrt(Ra*Pr).
-visc 0.026646
-diff_t 0.037538
-buoyancy true
-T_init 0.0

# Ra = 1e5, Pr = 0.71
# -visc 0.0026645825
# -diff_t 0.0037529465
# -n_timestep 4000
############################################################
############################################################

# Zone 4: left wall, hot (no-slip, u=v=0 by default)
-bc_4_type wall
-bc_4_T 1.0

# Zone 2: right wall, cold
-bc_2_type wall
-bc_2_T 0.0

# Zones 3, 5: bottom/top walls, adiabatic (no -bc_<zone>_T -> zero-gradient)
-bc_3_type wall
-bc_5_type wall

# Zone 7: interior faces (shared by two cells, no B.C.)
-bc_7_type interior

# No zone provides a pressure B.C. here, so the solver auto-anchors cell 0's pressure
# to zero (see PRESSURE_ANCHOR_NEEDED in inc_problems.C) to make the system solvable.
