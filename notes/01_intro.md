# 1. Introduction

Incompressible, single-phase flow is the working assumption behind a large share of
thermal-hydraulics (T/H) analysis for many nuclear systems, such as liquid metal cooled, molten salt cooled/fueled reactors.
For gas cooled reactors operating at low flow speed, the incompressible assumption is still appropriate.
Within that approximation, the physics of interest such as buoyancy-driven circulation, mixing, conjugate heat
transfer can be assumed to be governed entirely by the incompressible Navier-Stokes
equations coupled to an energy transport equation. Solving that system
accurately and robustly, on the unstructured meshes needed for realistic reactor geometries, is
of interest to the reactor thermal-hydraulics community.

The cell-centered, collocated finite-volume method (FVM) is the natural discretization choice
for this class of problem: it conserves mass, momentum, and energy exactly at the discrete level
on arbitrary unstructured meshes.
Its well-known complication is pressure-velocity coupling on a collocated mesh,
which can be addressed by the Rhie-Chow momentum interpolation.
This report works with a cell-centered FVM scheme, following the segregated predictor-corrector (SIMPLE-type)
algorithm of Dalal, Eswaran, and Biswas (DEB08).

Segregated algorithms have been the default in general-purpose CFD for decades.
They solve the momentum and pressure (and, where present, energy) equations in sequence
with proper linearizations and time lagging.
Compared with Newton's method, these algorithms are inexpensive and matrix-light.
On the other hand, these methods can be seen as
an approximate linearization: the segregated sweep can be shown to be Newton's
method applied to the fully coupled residual with part of the Jacobian dropped. An alternative is
to retain the full Jacobian and solve the coupled, nonlinear system monolithically with Newton's
method, using a Krylov linear solver (Newton-Krylov) for each Newton step. This removes the
lagging between fields entirely, at the cost of forming and preconditioning a much larger, fully
coupled linear system at every step.

For T/H applications, many of the flows of interest are the ones where segregated lagging could be problematic,
e.g., buoyancy-driven natural convection, in which momentum and energy are two-way coupled, is the
dominant mechanism in passive safety systems, and its transient response is exactly what a lagged,
one-step-behind coupling could be inaccurate.
A monolithic, fully implicit solver is also the natural building block for coupling a multidimensional CFD-like model
to the one-dimensional system-level models used for whole-plant safety analysis.
Tightly (implicitly) coupling a 1-D and a multi-D domain at a shared interface is far
more straightforward when the multi-D side already exposes a single coupled nonlinear residual and Jacobian.
Understanding what the monolithic formulation costs, and what it buys, in a standalone multidimensional solver
is a necessary step before using it that way.

We organize the remaining report with the following sections.
[Section 2](02_fvm.md) recaps the segregated DEB08 algorithm and derives the corresponding fully implicit,
monolithic Newton-Krylov formulation of the same discrete equations.
[Section 3](03_benchmarks.md) verifies both against standard incompressible-flow benchmarks,
such as lid-driven problem and flow past a cylinder problem, etc.
[Section 3.5](04_comparison.md) then compares the two methods in terms of numerical results and computational costs.
[Section 4](05_conclusion.md) provides a summary to this work and additional discussions to future work.
