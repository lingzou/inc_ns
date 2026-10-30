# 3.5 Segregated Solver vs. the Newton Solver

Section [2](02_fvm.md) showed that the segregated and monolithic formulations
are not identical.
This section provides numerical evidence to demonstrates that difference.
Section [3.5.2](#352-computational-cost) then shows that, despite the extra
cost of solving the coupled system, the Newton solver is practical in
practice: paired with a Schur-complement preconditioner it keeps the
slowdown relative to the segregated algorithm to a modest factor across
the benchmark suite.

## 3.5.1 The two solvers are not identical

In the segregated algorithm the mass flux is built from the
*provisional* velocity from the modified momentum equation (provisional momentum equation)
without the pressure gradient, and only after the flux-correction loop
has converged does the final re-solve (momentum equation) produce $u^{n+1}$.
The converged segregated flux is therefore $F_f = \overline{(u^{\ast})}_f\cdot\Delta S_f - \Delta t(\nabla p)_f\cdot\Delta S_f$,
while the monolithic flux (Rhie-Chow momentum interpolation, section [2.5](02_fvm.md#25-fully-implicit-monolithic-formulation)) is the same
expression written in terms of $u^{n+1}$, the velocity that actually solves the
momentum equation.

The two agree only to the extent that
$M^{-1}\Delta V\approx\Delta t\,I$, i.e.\ only to the extent that the
convection and diffusion parts of the momentum operator $M$ are
negligible against $\rho\Delta V/\Delta t$, i.e. very small time step size.
We demonstrate this using the lid-driven cavity problem with Re = 100,
and a coarse mesh with 830 cells.
Both solvers run to $t=0.2$,
the $L_2$ difference between the two final fields is

| $\Delta t$ | $u$ | $v$ | $p$ |
|:---:|:---:|:---:|:---:|
| 0.01   | $3.632\times10^{-2}$ | $3.918\times10^{-2}$ | $5.099\times10^{-2}$ |
| 0.005  | $1.224\times10^{-2}$ | $1.407\times10^{-2}$ | $2.223\times10^{-2}$ |
| 0.0025 | $4.311\times10^{-3}$ | $5.049\times10^{-3}$ | $1.052\times10^{-2}$ |

which clearly demonstrates that the two methods converge at smaller time step size.

It is noted that the above table is a whole-field norm.
Examining the data more carefully, it can be found that the differences concentrate in the top
corners, where the moving lid's prescribed velocity meets the stationary side walls.
This is why the segregated-versus-monolithic centerline comparison of section [3](03_benchmarks.md)
(Figures 3.4 and 3.5), sampled along $x=0.5$ and $y=0.5$ did not reveal the difference.

A second check was performed using the flow-past-a-cylinder benchmark of section [3](03_benchmarks.md).
The non-trivial differences in drags, a surface integrated quantity,
listed in Table 3.2, confirm that the two formulations are not identical.

## 3.5.2 Computational cost

The monolithic solve is more expensive per time step, and how much more depends on
how the Newton system is linearized and preconditioned. For each of the four benchmark problems of
section [3](03_benchmarks.md), run with the segregated algorithm and with four variants of the Newton-Krylov
solve of the Jacobian system (Jacobian structure in section [2.5.1](02_fvm.md#251-newton-linearization)) - an assembled or
matrix-free (mf) Jacobian action, each paired with either a Schur-complement or a plain direct (LU) preconditioner:

| case | segregated (s/step) | Schur (s/step) | mf+Schur (s/step) | LU (s/step) | mf+LU (s/step) |
|:---|:---:|:---:|:---:|:---:|:---:|
| cavity                | 0.00518 | 0.0536 | 0.0376 | 0.0676 | 0.0671 |
| heated/cooled cavity  | 0.00339 | 0.0727 | 0.0595 | 0.0981 | 0.0996 |
| cylinder              | 0.0607  | 2.296  | 0.893  | 1.709  | 1.714  |
| backward-facing step  | 0.204   | 2.146  | 2.944  | 3.853  | 3.771  |

From the table, it can be generally observed that the direct LU solver is the most expensive choice (as expected),
while the mf + Schur is the cheapest, with one exception that the flow past cylinder problem Schur is the most expensive.
For all Newton solver, Newton non-linear iteration counts
range from $2.0$ to $4.8$ per time step across every case and variant tested, essentially insensitive
to mesh size, Reynolds number, or the assembled/matrix-free choice (the LU and mf+LU variants, which
share the same exact preconditioner, converge in identical iteration counts in every case).

Where the cost actually goes depends on the preconditioner. For the Schur-preconditioned
variants (assembled or matrix-free), the finite-difference, coloring-based Jacobian evaluation used
to build the preconditioner dominates: $26$--$71\%$ of total solve time, with the Krylov solve taking
most of the remainder.
This implies that the Schur-preconditioned method could be further imporoved replacing the
expensive Jacobian evaluation with a hand-coded Jacobian evaluation subroutine.
For the LU-preconditioned variants, the dominant cost shifts instead
to the direct factorization of the full coupled Jacobian
$42$--$84\%$ of total solve time, and this fraction grows with problem size, from $42$--$44\%$ on the
${\sim}830$-cell cavity cases up to $69$--$84\%$ on the larger cylinder, backward-facing-step,
since direct factorization cost scales worse with problem size than the coloring-based Jacobian evaluation does;
correspondingly, Jacobian assembly itself falls to a minor $15$--$29\%$ of the total on those larger LU runs.

Even with the modest sizes of the problems, the Schur complement algorithm already demonstrate
superb performance than the LU direct solver,
which can only be even worser with larger problem size.
As discussed, the Schur-preconditioned method could be further imporoved using hand-coded Jacobian
evaluation subroutines,
indicating that the Schur-preconditioned method would be the more favored choice for larger-size problems.
