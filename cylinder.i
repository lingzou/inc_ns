-problem cylinder
-mesh flow-past-a-cylinder.msh
-dt 0.05
-visc 0.001
-n_timestep 20

# Zone 1: velocity inlet
-bc_1_type velocity_inlet
-bc_1_u 1.0
-bc_1_v 0.0

# Zones 2 and 4: horizontal (top/bottom) slip walls, v = 0
-bc_2_type horizontal_slip_wall
-bc_4_type horizontal_slip_wall

# Zone 3: pressure outlet
-bc_3_type pressure_outlet
-bc_3_p 0.0

# Zone 5: cylinder wall, no-slip
-bc_5_type wall

# Zone 7: interior faces (shared by two cells, no B.C.)
-bc_7_type interior
