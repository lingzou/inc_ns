# 4. Summary

In this work, we built, verified, and compared two ways of solving the incompressible,
cell-centered finite-volume equations on unstructured meshes, i.e., the segregated predictor-corrector
algorithm of Dalal, Eswaran, and Biswas (DEB08), and a fully implicit, monolithic Newton-Krylov alternative
extended from the DEB08 method.
In this work, the following has been achieved:

- implemented the DEB08 segregated (Picard/SIMPLE-type) algorithm with the optional
  Boussinesq-coupled temperature equation;
- derived the fully implicit monolithic formulation,
  extending DEB08's face-flux treatment to the Rhie-Chow style face flux
  appropriate for the monolithic formulation,
  and the pressure Schur-complement structure for preconditioning;
- implemented the monolithic formulation as a Newton-Krylov solve, with both assembled and
  matrix-free Jacobian actions and both Schur-complement and direct (LU) preconditioning;
- verified both solvers against standard incompressible-flow benchmarks,
  and against each other, confirming that both formulations are distinct
  but close and consistent discretizations that converge to each other at small time step size; and
- quantified monolithic solver costs compared with the segregated solver - roughly 10 times
  more expensive if picking the best-performance Newton variant, in general the mf + Schur algorithm,
  which could potentially be further improved using hand-coded Jacobian
  to drive down the cost associated with finite-differencing Jacobian formulations.

