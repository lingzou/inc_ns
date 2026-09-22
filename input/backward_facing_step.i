-problem backward_facing_step
-mesh input/backward_facing_step.msh
-dt 0.05
# See page 478 of Armaly et al., J. Fluid Mech., vol. 127, pp. 473–496, 1983.
# Re = u * (2*h) / nu, with u = avg inlet velocity, h = inlet channel height, nu = viscosity
# u = 1.5 [m/s], h = 1 [m], Re = 300 -> nu = 1.5 * 2 / 300 = 0.01
-visc 0.01
-n_timestep 20

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
# Average velocity of a parabolic profile is (2/3)*u_max, so u_max = 1.5 * 1.5 = 2.25
# keeps the same avg inlet velocity (and thus Re = 300) as the uniform profile above.
-bc_1_type velocity_inlet
-bc_1_u_profile parabolic
-bc_1_u_ymin -1.0
-bc_1_u_ymax 0.0
-bc_1_u_umax 2.25
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
