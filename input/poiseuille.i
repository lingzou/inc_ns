-problem poiseuille
-mesh input/poiseuille.msh
# Plane Poiseuille flow: channel of length L = 10 [m], plate spacing h = 1 [m],
# uniform inlet u = U = 1 [m/s], outlet p = 0, no-slip plates.
# Re = U * h / nu = 1 * 1 / 0.01 = 100
-visc 0.01
-dt 0.1
-n_timestep 50
-output_interval 10

# Zone 1: inlet -- uniform velocity inlet
-bc_1_type velocity_inlet
-bc_1_u 1.0
-bc_1_v 0.0

# Zones 2, 4: bottom_wall, top_wall -- the two plates, solid no-slip
-bc_2_type wall
-bc_4_type wall

# Zone 3: outlet -- pressure outlet
-bc_3_type pressure_outlet
-bc_3_p 0.0

# Zone 7: interior faces (shared by two cells, no B.C.)
-bc_7_type interior
