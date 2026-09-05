-problem cavity
-mesh input/cavity.msh
-dt 0.01
-visc 0.01
-n_timestep 20

# Zone 5: lid -- wall with prescribed (moving) velocity
-bc_5_type wall
-bc_5_u 1.0
-bc_5_v 0.0

# Zones 2, 3, 4: stationary walls -- -bc_<zone>_u/_v omitted, default to 0
-bc_2_type wall
-bc_3_type wall
-bc_4_type wall

# Zone 7: interior faces (shared by two cells, no B.C.)
-bc_7_type interior

# No zone provides a pressure B.C. here, so the solver auto-anchors cell 0's pressure
# to zero (see PRESSURE_ANCHOR_NEEDED in inc_problems.C) to make the system solvable.
