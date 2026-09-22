-problem cylinder
-mesh input/flow-past-a-cylinder.msh
-dt 0.05
-visc 0.001
-n_timestep 20

# DEB08's convective blend factor: 1.0 = pure central difference, 0.0 = pure first-order
# upwind. Set explicitly rather than inheriting the 0.5 default -- 50% upwinding is a lot
# of numerical diffusion to apply silently, and on the separated flows it matters (it
# shortens the backward-facing step's recirculation bubble by 16% at Re = 300). At Re = 40
# it also moves this case's drag: C_D = 1.810 at gamma = 0.5 vs 1.729 at gamma = 1.0,
# against DEB08 Table 2's 1.599.
-gamma 1.0

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
