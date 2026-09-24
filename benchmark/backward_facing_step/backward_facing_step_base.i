-problem backward_facing_step
-mesh input/backward_facing_step.msh
-dt 0.05
# See page 478 of Armaly et al., J. Fluid Mech., vol. 127, pp. 473–496, 1983.
# Re = u_avg * (2*h) / nu, with u_avg = avg inlet velocity, h = inlet channel height (both
# 1.0 here, see the inlet B.C. below) -> nu = 2 / Re. -visc alone sets Re; see
# run_Re_sweep.py, which sweeps Re = 50, 100, 150, 200, 250, 300, 400 over this file.
# This default (nu = 0.01) is the Re = 200 case.
-visc 0.01
# Same dt / n_timestep for every Re in the sweep. Only Re = 200 (this default) and
# Re = 300 (input/backward_facing_step.i's u_max = 2.25 case) have been checked against
# Armaly et al. -- inspect each case's "<name>_out_convergence.csv" (u_L2_error column)
# before trusting reattachment lengths at the other Re values.
-n_timestep 2000
-output_interval 100

# DEB08's convective blend factor: 1.0 = pure central difference, 0.0 = pure first-order
# upwind. Set explicitly rather than inheriting the 0.5 default. This case is the most
# gamma-sensitive one in the repo: run to steady state (t = 150), the reattachment length
# is Lr/S = 5.64 at gamma = 0.5 but 6.65 at gamma = 1.0, against Armaly et al.'s ~6.6 --
# so the default alone accounts for a 16% under-prediction. Reattachment is set by the
# spreading rate of the free shear layer off the step, which false diffusion attacks
# directly; see docs/newton_krylov_formulation.tex section "Choice of the blend factor".
-gamma 1.0

# Zone 1: inlet -- velocity inlet, parabolic profile
# Inlet spans y in [-1, 0] (see mesh/backward_facing_step.geo, points A/B).
# Average velocity of a parabolic profile is (2/3)*u_max, so u_max = 1.5 gives u_avg = 1.0.
# This is fixed across the Re sweep -- -visc alone sets Re, not the inlet velocity.
-bc_1_type velocity_inlet
-bc_1_u_profile parabolic
-bc_1_u_ymin -1.0
-bc_1_u_ymax 0.0
-bc_1_u_umax 1.5
-bc_1_v 0.0

# Zones 2, 3, 5: step, bottom_wall, top_wall -- solid no-slip walls
-bc_2_type wall
-bc_3_type wall
-bc_5_type wall

# Zone 4: outlet -- pressure outlet
-bc_4_type pressure_outlet
-bc_4_p 0.0

# Zone 7: interior faces (shared by two cells, no B.C.)
-bc_7_type interior
