# Notes on Cell-centered Finite Volume Method and Newton-Krylov Method in Solving Incompressible Flow Problems

**Author:** Ling Zou
**Date:** Oct. 30, 2026

## Abstract

Incompressible single-phase flow is fundamental to thermal-hydraulics analysis encountered in many nuclear systems, such as liquid metal cooled, molten salt cooled/fueled reactors.
This report develops, implements, and compares two cell-centered, collocated finite-volume solvers for this class of problem on unstructured meshes: the segregated predictor-corrector algorithm of Dalal, Eswaran, and Biswas (DEB08), and a fully implicit, monolithic Newton-Krylov alternative derived from the same scheme.
In the monolithic formulation, we extend the DEB08's face flux treatment using the Rhie-Chow momentum interpolation.
In addition, a pressure Schur complement method was derived leveraging the segregated solver knowledge.
Both solvers are implemented with the Boussinesq-coupled energy equation and verified against standard incompressible-flow benchmarks, including the lid-driven cavity, backward-facing step, flow past a cylinder, and the differentially heated cavity problems.
The two formulations are then compared and it is clear that they are distinct consistent discretizations that converge to each other under time-step refinement.
Numerical tests showed that the Newton's method converges reliably for all problems but more expensive (as expected).
The Schur complement method has also demonstrated great potential compared with the direct LU solver option.

## Contents

1. [Introduction](01_intro.md)
2. [Dalal-Eswaran-Biswas (DEB08) Method](02_fvm.md)
   - 2.1 Governing Equations
   - 2.2 Discretized Form
     - 2.2.1 Convection term
     - 2.2.2 Gradient term
     - 2.2.3 Diffusion term
   - 2.3 Segregated Solution Algorithm
   - 2.4 Additional Notes
   - 2.5 Fully Implicit (monolithic) Formulation
     - 2.5.1 Newton linearization
     - 2.5.2 Pressure Schur complement
   - 2.6 Including a Temperature Field
3. [Code Benchmarks](03_benchmarks.md)
   - 3.1 Lid Driven Problem
   - 3.2 Backward-facing Step Flow
   - 3.3 Flow Past a Cylinder
   - 3.4 Differentially Heated Cavity
   - 3.5 [Segregated Solver vs. the Newton Solver](04_comparison.md)
     - 3.5.1 The two solvers are not identical
     - 3.5.2 Computational cost
4. [Summary](05_conclusion.md)
   - [Acknowledgments](05_conclusion.md#acknowledgments)
5. [References](references.md)
