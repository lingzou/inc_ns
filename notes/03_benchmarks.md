# 3. Code Benchmarks

## 3.1 Lid Driven Problem

The lid-driven cavity is the one of the standard tests for an incompressible flow solver.
The domain is a unit square, $`[0,1]\times[0,1]`$.
The top wall moves at $`u=1`$, $`v=0`$; the other three walls are stationary and no-slip.
Pressure is pinned at a single reference cell since no boundary admits a Dirichlet pressure condition.
In this work, we performed transient simulations for $`Re=100`$, $`400`$, $`1000`$, and $`3200`$ to the steady state,
using the same unstructured triangular mesh of 8238 cells, as shown in Figure 3.1.
Each side contains roughly 64 cells.
Figure 3.2 shows the stream function of each case,
which compares well against the classical benchmark solutions of Ghia, et al. (Ghia et al., 1982):
the primary vortex center and the secondary corner eddies that appear at the higher Reynolds numbers
match the expected locations and relative sizes.

<p align="center">
  <img src="figs/lid_driven/lid_driven_mesh.png" alt="Lid-driven mesh" width="45%" />
</p>

*Figure 3.1: Lid driven problem: unstructured triangular mesh (8238 cells) obtained from a
$`64\times64`$ transfinite grid with alternating diagonals.*

<p align="center">
  <img src="figs/lid_driven/cavity_Re100_streamfunction.png" alt="Re = 100" width="45%" />
  <img src="figs/lid_driven/cavity_Re400_streamfunction.png" alt="Re = 400" width="45%" />
  <br/>
  <img src="figs/lid_driven/cavity_Re1000_streamfunction.png" alt="Re = 1000" width="45%" />
  <img src="figs/lid_driven/cavity_Re3200_streamfunction.png" alt="Re = 3200" width="45%" />
</p>

*Figure 3.2: Lid driven problem: stream function contours at Re = 100, 400, 1000, and 3200.*

A more quantitative comparison is obtained by extracting the velocity profiles along the two
centerlines of the cavity and comparing them directly against the tabulated data of Ghia, et al. (Ghia et al., 1982),
as shown in Figure 3.3.
The $`u`$-velocity along the vertical centerline $`x=0.5`$ and the $`v`$-velocity along the horizontal centerline $`y=0.5`$
both show excellent agreement with the reference data at all four Reynolds numbers.
As $`Re`$ increases, the profiles develop the characteristic S-shape with steeper near-wall gradients and
larger-magnitude peaks, reflecting the thinning boundary layers and the strengthening secondary
eddies. The agreement is closest at $`Re=100`$ and $`400`$, where the boundary layers remain well
resolved by the mesh; at $`Re=1000`$ and especially $`Re=3200`$ small but visible deviations appear near
the walls where the velocity gradients are steepest.
Mesh refinement study confirmed that they can be resolved better using finer meshes,
which are not included in this report.

<p align="center">
  <img src="figs/lid_driven/u_mid_comparison.png" alt="u along x = 0.5" width="48%" />
  <img src="figs/lid_driven/v_mid_comparison.png" alt="v along y = 0.5" width="48%" />
</p>

*Figure 3.3: Lid driven problem: centerline velocity profiles compared against Ghia, Ghia, and
Shin (Ghia et al., 1982) for Re = 100, 400, 1000, and 3200. (a) $`u`$ along the vertical centerline ($`x=0.5`$);
(b) $`v`$ along the horizontal centerline ($`y=0.5`$).*

The same centerline profiles also give a direct check of the segregated solver against the
monolithic Newton-Krylov solver (section [2.5](02_fvm.md#25-fully-implicit-monolithic-formulation)) on this mesh, both at steady state and
at time $`t=0.2\,\mathrm{s}`$ (20 steps at $`\Delta t=0.01`$ from t = 0).
On the 8238-cell mesh used throughout this section, Figures 3.4 and 3.5 show the two solvers' profiles
essentially on top of each other at all four Reynolds numbers, both at steady state and at $`t=0.2\,\mathrm{s}`$.
This visual agreement should not be read as the two solvers being identical.
As discussed in section [3.5.1](04_comparison.md#351-the-two-solvers-are-not-identical), they do not solve exactly
the same discrete equations: the segregated flux is built from the momentum predictor rather than
the fully re-solved velocity. A whole-field comparison on this same case (not just the two centerlines
plotted here) does show a measurable difference, concentrated mainly at the two top corners where the
moving lid meets the stationary side walls, see more details in section [3.5.1](04_comparison.md#351-the-two-solvers-are-not-identical).

<p align="center">
  <img src="figs/lid_driven/u_mid_solver_comparison.png" alt="u along x = 0.5, steady" width="48%" />
  <img src="figs/lid_driven/v_mid_solver_comparison.png" alt="v along y = 0.5, steady" width="48%" />
</p>

*Figure 3.4: Lid driven problem: segregated versus monolithic (Newton) centerline velocity profiles
at steady state, for Re = 100, 400, 1000, and 3200. (a) $`u`$ along the vertical centerline ($`x=0.5`$);
(b) $`v`$ along the horizontal centerline ($`y=0.5`$).*

<p align="center">
  <img src="figs/lid_driven/u_mid_t0.2_solver_comparison.png" alt="u along x = 0.5, t = 0.2 s" width="48%" />
  <img src="figs/lid_driven/v_mid_t0.2_solver_comparison.png" alt="v along y = 0.5, t = 0.2 s" width="48%" />
</p>

*Figure 3.5: Lid driven problem: segregated versus monolithic (Newton) centerline velocity profiles
at $`t=0.2\,\mathrm{s}`$, for Re = 100, 400, 1000, and 3200. (a) $`u`$ along the vertical centerline ($`x=0.5`$);
(b) $`v`$ along the horizontal centerline ($`y=0.5`$).*

## 3.2 Backward-facing Step Flow (DEB08)

The backward-facing step adds a geometric singularity, the step corner, that produces a
separated shear layer and a recirculation zone whose reattachment point is a sensitive,
well-documented measure of solver accuracy. The domain consists of an inlet duct of length = 1 and height $`h=1`$
($`0\le x\le 1`$, $`-1\le y\le 0`$) that expands abruptly at $`x=1`$ into a channel of height $`h+S`$,
with the step of height $`S=0.9423`$ taking the outlet-to-inlet height ratio to $`1.9423`$,
matching the expansion ratio used in the experiments of Armaly et al. (Armaly et al., 1983); the
downstream channel extends to $`x=21`$ to let the recirculation and reattachment develop fully.
A fully-developed parabolic velocity profile with a fixed bulk average $`u_{\mathrm{avg}}=1`$ is imposed at the inlet,
all walls are no-slip, and the outlet is a zero-gauge-pressure boundary.
The Reynolds number is defined on the inlet duct's hydraulic diameter, $`Re = u_{\mathrm{avg}}(2h)/\nu`$.
Figure 3.6 show the velocity results for different Re.
The mesh has 21076 triangular cells, intentionally generated fine near the step and coars downstream.
Because upwind-biased convection smears the thin separated shear layer and systematically
shortens the recirculation bubble, these cases use the central-differencing blend ($`\gamma=1.0`$).
With that setting the predicted reattachment length $`L_r`$, normalized by the step height $`S`$, agrees well
reference values (see discussion below).

<p align="center">
  <img src="figs/backward_facing_step/mesh.png" alt="Backward-facing step mesh" width="90%" />
  <br/>
  <img src="figs/backward_facing_step/Re100.png" alt="Re = 100" width="90%" />
  <br/>
  <img src="figs/backward_facing_step/Re200.png" alt="Re = 200" width="90%" />
  <br/>
  <img src="figs/backward_facing_step/Re300.png" alt="Re = 300" width="90%" />
  <br/>
  <img src="figs/backward_facing_step/Re400.png" alt="Re = 400" width="90%" />
</p>

*Figure 3.6: Backward-facing step flow: computational mesh and velocity fields at Re = 100, 200, 300, and 400.*

A more systematic check on reattachement length was then performed using two-dimensional
backward-facing-step reference study of Biswas et al. (Biswas et al., 2004) for the same expansion ratio, 1.9423.
Table 3.1 and Figure 3.7 compare this study with the Biswas reference results.
It can be found that the present solver matches Biswas et al.'s reattachment lengths very well at every Reynolds number.

**Table 3.1:** Backward-facing step flow: normalized reattachment length $`L_r/S`$ compared against
Biswas, Breuer, and Durst (Biswas et al., 2004).

| $`Re`$ | $`L_r/S`$ (Biswas et al.) | $`L_r/S`$ (this work) | Error (%) |
|:---:|:---:|:---:|:---:|
| 50  | 1.686 | 1.678 | $`-0.47`$ |
| 100 | 2.872 | 2.852 | $`-0.69`$ |
| 150 | 3.910 | 3.912 | $`0.04`$  |
| 200 | 4.889 | 4.877 | $`-0.25`$ |
| 250 | 5.824 | 5.780 | $`-0.76`$ |
| 300 | 6.664 | 6.646 | $`-0.27`$ |
| 400 | 8.172 | 8.152 | $`-0.24`$ |

<p align="center">
  <img src="figs/backward_facing_step/comparison.png" alt="Reattachment length vs Re" width="60%" />
</p>

*Figure 3.7: Backward-facing step flow: normalized reattachment length $`L_r/S`$ versus Reynolds number, compared against
Biswas, Breuer, and Durst (Biswas et al., 2004).*

## 3.3 Flow Past a Cylinder

Flow past a circular cylinder is an external-flow benchmark, in contrast to the two enclosed-domain
cases above, and tests the solver's handling of an immersed no-slip boundary together with a
far-field that must be large enough not to contaminate the computed loads. The domain is a full
rectangle $`50\times25`$ with a cylinder of diameter $`D=1`$ centered at $`(25,12.5)`$,
i.e., 25 diameters from both the inlet and the outlet, giving a blockage ratio of only 4\% so that
the top and bottom boundaries do not meaningfully confine the flow. A uniform velocity $`u=1`$ is
imposed at the inlet, the top and bottom boundaries are slip walls, the outlet is a zero-gauge-pressure
boundary, and the cylinder surface is no-slip, with the drag and lift forces obtained by integrating
the traction over that surface. The Reynolds number is $`Re = UD/\nu = 1/\nu`$, and $`Re=10`$, $`20`$, and
$`40`$ all fall in the steady, symmetric-wake regime, below the onset of
vortex shedding. The mesh has 24874 triangular cells, with 80 elements around the cylinder
circumference, a geometric growth ratio of 1.05 away from the surface, and a refined wake corridor
extending 14.5 diameters downstream. The velocity and pressure fields shown in
Figures 3.8 and 3.9, together with the integrated pressure ($`C_p`$) and skin-friction ($`C_f`$) contributions to the drag
coefficient $`C_D=C_p+C_f`$, are summarized in Table 3.2 alongside the reference values of Dennis and Chang
(Dennis \& Chang, 1970) and of DEB08 (DEB08).
The computed $`C_D`$ agrees with Dennis and Chang to within $`5.8\%`$ and with DEB08 to within
$`3.2\%`$ over this Reynolds number range.
As $`Re`$ increases from $`10`$ to $`40`$, the pressure fraction of the drag, $`C_p/C_D`$, rises from
$`0.553`$ to $`0.661`$, reflecting the thinning boundary layer and the correspondingly shrinking
relative contribution of skin friction; the computed lift coefficient $`C_l`$ stays at or below
$`0.002`$ at all three Reynolds numbers, consistent with the expected steady, symmetric wake below
the onset of vortex shedding.

**Table 3.2:** Flow past a cylinder: drag coefficient breakdown compared against Dennis and
Chang (Dennis \& Chang, 1970) and DEB08 (DEB08), and between the two solvers of this work.

| | | $`C_f`$ | $`C_p`$ | $`C_D`$ |
|:---|:---|:---:|:---:|:---:|
| **$`Re = 10`$** | This work (segregated) | 1.284 | 1.589 | 2.874 |
| | This work (Newton)     | 1.294 | 1.630 | 2.925 |
| | DEB08                  | 1.288 | 1.680 | 2.968 |
| | Dennis \& Chang        | 1.246 | 1.600 | 2.846 |
| **$`Re = 20`$** | This work (segregated) | 0.838 | 1.275 | 2.113 |
| | This work (Newton)     | 0.847 | 1.270 | 2.118 |
| | DEB08                  | 0.838 | 1.308 | 2.146 |
| | Dennis \& Chang        | 0.812 | 1.233 | 2.045 |
| **$`Re = 40`$** | This work (segregated) | 0.546 | 1.064 | 1.611 |
| | This work (Newton)     | 0.555 | 1.043 | 1.598 |
| | DEB08                  | 0.537 | 1.062 | 1.599 |
| | Dennis \& Chang        | 0.522 | 0.998 | 1.522 |

Table 3.2 also runs the same three cases through the monolithic Newton-Krylov solver (matrix-free, Schur-complement
preconditioner), each driven to the same steady state as its segregated counterpart. The two solvers' $`C_D`$ values are close but
not identical. The differences are far larger than either solver's own residual at convergence, so they are not
solver-tolerance noise; they are the same segregated/monolithic discrepancy quantified directly in
section [3.5.1](04_comparison.md#351-the-two-solvers-are-not-identical), here showing up as a difference in an integrated
whole-surface engineering quantity.

<p align="center">
  <img src="figs/flow_past_cylinder/mesh_around_cylinder.png" alt="Cylinder mesh" width="90%" />
  <br/>
  <img src="figs/flow_past_cylinder/vel_Re10.png" alt="Velocity, Re = 10" width="90%" />
  <br/>
  <img src="figs/flow_past_cylinder/vel_Re20.png" alt="Velocity, Re = 20" width="90%" />
  <br/>
  <img src="figs/flow_past_cylinder/vel_Re40.png" alt="Velocity, Re = 40" width="90%" />
</p>

*Figure 3.8: Flow past a cylinder: computational mesh around the cylinder and velocity fields at Re = 10, 20, and 40.*

<p align="center">
  <img src="figs/flow_past_cylinder/cylinder_Re10_out_pressure.png" alt="Pressure, Re = 10" width="45%" />
  <img src="figs/flow_past_cylinder/cylinder_Re20_out_pressure.png" alt="Pressure, Re = 20" width="45%" />
  <br/>
  <img src="figs/flow_past_cylinder/cylinder_Re40_out_pressure.png" alt="Pressure, Re = 40" width="45%" />
</p>

*Figure 3.9: Flow past a cylinder: pressure fields at Re = 10, 20, and 40.*

## 3.4 Differentially Heated Cavity

This case uses the same unit-square domain and mesh as the lid-driven cavity, Figure 3.1, but removes the
moving lid and instead drives the flow by buoyancy, exercising the Boussinesq source term and the
momentum-energy coupling that the two previous isothermal cases never touch. All four walls are
no-slip; the left wall is held hot at $`T=1`$, the right wall cold at $`T=0`$, and the top and bottom
walls are adiabatic (zero normal temperature gradient). The Prandtl number is fixed at $`Pr=0.71`$
(air), and the Rayleigh number is set by scaling the viscosity and thermal diffusivity together,
$`\nu = \sqrt{Pr/Ra}`$ and $`\alpha = 1/\sqrt{Ra\cdot Pr}`$ (see DEB08 (DEB08)).
In this work, we studided the cases for $`Ra = 10^3`$, $`10^4`$, $`10^5`$, and $`10^6`$.
Figure 3.10 shows the temperature fields of these cases, spanning conduction-dominated transport at the lowest
$`Ra`$ to a boundary-layer-dominated convective regime at the highest.
The steady-state average Nusselt numbers, evaluated separately on the hot
and cold walls as $`\overline{Nu}_{\text{hot}}`$ and $`\overline{Nu}_{\text{cold}}`$ (at an exactly
energy-conserving steady state the two should coincide)
are reported in Table 3.3.
The mismatch between the hot and code side values stays below $`0.05\%`$ at every $`Ra`$.
The computed $`\overline{Nu}`$ is compared against the classical benchmark of de Vahl
Davis (de Vahl Davis, 1983) and against DEB08 (DEB08).
Agreement with both references are very well, while it is observed that errors are largest
for the largest $`Ra`$ case, most likely due to the coarse mesh being used not able to resolve the thin
thermal boundary layers.

**Table 3.3:** Differentially heated cavity: average Nusselt number on the hot and cold walls,
their mismatch, and the resulting $`\overline{Nu}`$, compared against de Vahl
Davis (de Vahl Davis, 1983) and DEB08 (DEB08), alongside the $`L^2`$ norm of
the steady-state temperature residual.

| $`Ra`$ | $`\overline{Nu}_{\text{hot}}`$ | $`\overline{Nu}_{\text{cold}}`$ | Mismatch (%) | $`\overline{Nu}`$ | DEB08 | Error (%) | de Vahl Davis | Error (%) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $`10^3`$ | 1.1222 | 1.1219 | 0.032 | 1.1220 | 1.160 | $`-3.27`$ | 1.118 | 0.36 |
| $`10^4`$ | 2.2547 | 2.2535 | 0.050 | 2.2541 | 2.247 | $`0.32`$  | 2.243 | 0.49 |
| $`10^5`$ | 4.5938 | 4.5924 | 0.030 | 4.5931 | 4.613 | $`-0.43`$ | 4.519 | 1.64 |
| $`10^6`$ | 9.2848 | 9.2835 | 0.014 | 9.2842 | 9.106 | $`1.96`$  | 8.800 | 5.51 |

<p align="center">
  <img src="figs/heated_cavity/cavity_buoyancy_Ra1e3_temperature.png" alt="Ra = 1e3" width="45%" />
  <img src="figs/heated_cavity/cavity_buoyancy_Ra1e4_temperature.png" alt="Ra = 1e4" width="45%" />
  <br/>
  <img src="figs/heated_cavity/cavity_buoyancy_Ra1e5_temperature.png" alt="Ra = 1e5" width="45%" />
  <img src="figs/heated_cavity/cavity_buoyancy_Ra1e6_temperature.png" alt="Ra = 1e6" width="45%" />
</p>

*Figure 3.10: Differentially heated cavity problem: temperature contours at Ra = $`10^3`$, $`10^4`$, $`10^5`$, and $`10^6`$.*
