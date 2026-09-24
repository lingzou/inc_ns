-problem cylinder
-mesh benchmark/flow_pass_cylinder/flow-pass-a-cylinder.msh
-dt 0.05
-visc 0.1
-n_timestep 1400
-output_interval 100

# D = 1 and U = 1 (zone 1 below), so Re = U*D/nu = 1/nu and -visc above is the
# only Re knob: 0.1 here is Re = 10. run_Re_sweep.py overwrites -visc (and, where
# needed, -n_timestep/-output_interval) to generate the Re = 10/20/40 cases.
# All three are below the ~Re = 47 vortex-shedding threshold, so each wake stays
# steady; 1400 steps (t = 70) is enough to reach it at Re = 40, the slowest.

# DEB08's convective blend factor: 1.0 = pure central difference, 0.0 = pure first-order
# upwind. Set explicitly rather than inheriting the 0.5 default -- 50% upwinding is a lot
# of numerical diffusion to apply silently, and on the separated flows it matters (it
# shortens the backward-facing step's recirculation bubble by 16% at Re = 300). At Re = 40
# it also moves this case's drag: C_D = 1.810 at gamma = 0.5 vs 1.729 at gamma = 1.0,
# against DEB08 Table 2's 1.599 (a known, tracked +13% discrepancy --
# docs/to-do-list-09-22-2026.md, item 13).
-gamma 0.5

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

# Surface force monitors on the cylinder (zone 5): appends Cf, Cp, Cd, Cl
# columns to the convergence CSV, computed each time step. Reference
# velocity/direction are taken from zone 1's inlet B.C.; reference length
# is zone 5's bounding-box size (~1.0, the cylinder diameter). Not
# compatible with -newton true.
-force_zones 5
