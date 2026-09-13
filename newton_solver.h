#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include <fstream>
#include <vector>

#include <petscsnes.h>

#include "FluentTwoDMesh.h"
#include "ExodusOutput.h"
#include "inc_problems.h"

/*
 * Monolithic Newton--Krylov time loop, selected by "-newton true" (see NEWTON_ON in
 * inc_problems.h). This is an alternative to -- not a replacement for -- the segregated
 * Picard/SIMPLE loop in cxx_main.C, which remains the default.
 *
 * Where the segregated path drives the face mass flux F_f to a fixed point with up to
 * 100 sub-iterations and then solves energy once, outside that loop, this path forms a
 * single nonlinear residual R(X) = 0 over the interlaced state X = (u,v,p,T) per cell
 * and hands it to PETSc's SNES. The derivation, the 4x4 Jacobian block structure, and
 * the reason this matters (the segregated energy lag makes -buoyancy true only a
 * pseudo-transient) are in docs/newton_krylov_formulation.tex.
 *
 * The residual deliberately reuses the operators cxx_main.C/inc_problems.C already
 * assemble rather than re-deriving the discretization: diffusion (including the
 * nonorthogonal cross terms) and the time-derivative diagonal are linear and already
 * live in M_*_FIXED, and updateAdvectionOperator() already builds the convective
 * operator from a given F_f. So, per cell,
 *
 *     R_u = M_USTAR(F_f(X)) * u - b_u(F_f(X))
 *     R_v = M_VSTAR(F_f(X)) * v - b_v(F_f(X), T)      // b_v carries the buoyancy source
 *     R_T = M_T(F_f(X))     * T - b_T(F_f(X))
 *     R_p = (1/dt) * sum_f F_f(X)                      // discrete continuity
 *
 * with F_f evaluated at the current iterate instead of a lagged one. Because both paths
 * then go through the same operators, with -buoyancy false they solve provably the same
 * discrete problem and must agree to solver tolerance -- see runNewtonTests.sh.
 */

struct NewtonSolverConfig
{
  FluentTwoDMesh * p_mesh;

  // Green-Gauss pressure-gradient reconstruction coefficients, precomputed once in
  // cxx_main.C. Only grad_p is needed here: the u/v/T gradient stencils exist solely to
  // build the nonorthogonal cross-diffusion into M_*_FIXED at assembly time, and that is
  // already baked in by the time this runs.
  std::vector<GRAD> * grad_p;

  // Assembled once in cxx_main.C, borrowed read-only: viscous (or thermal) diffusion
  // including nonorthogonal cross terms, plus the rho*V/dt time-derivative diagonal.
  Mat M_USTAR_FIXED, M_VSTAR_FIXED, M_T_FIXED;
  // Working copies with the same nonzero pattern. Rebuilt from the *_FIXED matrices on
  // every residual evaluation, so their contents on entry and exit are meaningless.
  Mat M_USTAR, M_VSTAR, M_T;

  // cxx_main.C's persistent RHS baselines, which at this point hold
  // (Dirichlet-BC diffusion contributions) + phi^0 * rho*V/dt. runNewtonSolver()
  // subtracts the phi^0 part back out once at entry to recover the genuinely
  // time-independent BC part, then rebuilds b = b_bc + phi^n * rho*V/dt itself at the
  // top of each step. That keeps phi^n structurally frozen across the Newton iterations
  // within a step, rather than relying on the in-place roll-forward the segregated loop
  // uses (cxx_main.C, end of the time loop).
  Vec b_USTAR, b_VSTAR, b_T;

  // Solution state, owned by cxx_main.C. On entry these hold the initial condition; on
  // return, the final time level, so cxx_main.C's cleanup and any post-processing work
  // unchanged.
  Vec u, v, p, T;
  // Face mass flux, sized n_Faces. Filled from the converged state each step so the
  // final-step stream function in writeExodusStep() still has what it needs.
  Vec F_face;

  int n_timestep;
  int output_interval;
  ExodusOutputContext * exoCtx;
  std::ofstream * convergence_file;
};

// Runs the whole time loop. Output (Exodus records, the convergence CSV, the
// output_interval / contiguous-record rule, the final-step stream function) mirrors the
// segregated loop's exactly, so the two modes' Exodus files are directly comparable.
void runNewtonSolver(NewtonSolverConfig & cfg);

#endif // NEWTON_SOLVER_H
