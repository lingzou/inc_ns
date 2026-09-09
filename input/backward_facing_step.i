-problem backward_facing_step
-mesh input/backward_facing_step.msh
-dt 0.05
# See page 478 of Armaly et al., J. Fluid Mech., vol. 127, pp. 473–496, 1983.
# Re = u * (2*h) / nu, with u = avg inlet velocity, h = inlet channel height, nu = viscosity
# u = 1.5 [m/s], h = 1 [m], Re = 300 -> nu = 1.5 * 2 / 300 = 0.01
-visc 0.01
-n_timestep 20

# Zone 1: inlet -- velocity inlet
-bc_1_type velocity_inlet
-bc_1_u 1.5
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
