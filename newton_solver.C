#include <cmath>
#include <iostream>
#include <set>

#include "newton_solver.h"

/*
 * See newton_solver.h for what this is and why the residual is built out of the existing
 * operators rather than re-deriving the discretization.
 *
 * All four blocks (u, v, p, T) are live. The energy equation is inside the nonlinear
 * solve rather than lagged behind it, which is the point of the exercise: it closes the
 * T -> buoyancy -> momentum -> advection -> T cycle that the segregated path leaves open.
 *
 * One thing NOT to expect: agreement with the segregated solver at fixed dt. The two are
 * different discrete problems -- the segregated scheme interpolates the face mass flux
 * from the predictor velocity u* (DEB-9 omits the pressure gradient), this one from the
 * velocity that actually solves momentum -- and they agree only to O(dt). See the
 * "...but they do not share a fixed point" subsection of
 * docs/newton_krylov_formulation.tex, and runNewtonTests.sh, which checks the
 * refinement rate rather than an equality.
 */

// DOF layout: interlaced, 4 per cell -- (u,v,p,T) for cell 0, then cell 1, ...
// Interlacing keeps each cell's 4x4 coupling block contiguous, which is what ILU and a
// later PCFIELDSPLIT both want. Cell ids from the mesh are 1-based; DOF indices are not.
static const PetscInt N_DOF_PER_CELL = 4;
static inline PetscInt dof_u(long int cell_id) { return N_DOF_PER_CELL * (cell_id - 1) + 0; }
static inline PetscInt dof_v(long int cell_id) { return N_DOF_PER_CELL * (cell_id - 1) + 1; }
static inline PetscInt dof_p(long int cell_id) { return N_DOF_PER_CELL * (cell_id - 1) + 2; }
static inline PetscInt dof_T(long int cell_id) { return N_DOF_PER_CELL * (cell_id - 1) + 3; }

namespace
{

// Everything FormResidual() needs. Pointers into cfg are borrowed, not owned; the Vecs
// and Mats created in runNewtonSolver() are destroyed there.
struct NewtonCtx
{
  NewtonSolverConfig * cfg;

  // Per-field scratch, refilled from X on every residual evaluation.
  Vec u_w, v_w, p_w, T_w;
  Vec grad_p_x, grad_p_y;
  Vec F_0f, F_f;              // face fluxes, sized n_Faces
  Vec b_u_work, b_v_work, b_T_work;

  // Right-hand-side baseline for the CURRENT time step:
  //   b_*_step = (Dirichlet-BC diffusion contributions) + phi^n * rho*V/dt
  // Rebuilt once per time step, held fixed across every Newton iteration within it.
  Vec b_u_step, b_v_step, b_T_step;
  Vec b_u_bc, b_v_bc, b_T_bc;  // the time-independent B.C. part alone

  // Old time level, frozen for the duration of a step.
  Vec u_old, v_old, T_old;

  // Residual scratch. Persistent rather than VecDuplicate'd per call: FD colouring
  // calls FormResidual() once per colour per Newton iteration, so an allocation in
  // there would be thousands of malloc/free pairs per time step.
  Vec R_u, R_v, R_T;
  std::vector<PetscScalar> div;
};

// Combined stand-in for updateMassVeclocities() + updateFfaceStar(), evaluated at the
// current iterate. We deliberately do NOT call updateMassVeclocities() itself: it also
// accumulates the divergence into a pressure RHS and force-zeroes bb_p[0] when the
// pressure is anchored, and that hidden side effect has no place inside a residual.
//
// The interpolation below uses face->distance_ratio(), matching updateMassVeclocities().
// That looks like a different weighting from the DEB-5 volume weighting alpha_12 used
// for the convected face value, but on a triangular mesh the two are the SAME NUMBER:
// the centroid-to-edge distance of a triangle is d = 2*volume/(3*edge_length), and the
// two cells share the edge, so d1/d2 = v1/v2 and 1 - d1/(d1+d2) = v2/(v1+v2) = alpha_12.
// Verified to roundoff (~1e-12) on cavity, cylinder and backward_facing_step. The
// identity is specific to triangles -- it would break if the mesh ever gained quads.
void computeFaceFluxes(FluentTwoDMesh * p_mesh, Vec u, Vec v, Vec p,
                       Vec gradP_x, Vec gradP_y, Vec F_0f, Vec F_f)
{
  const PetscScalar * uu, * vv, * pp, * gpx, * gpy;
  PetscScalar * f0f, * ff;
  VecGetArrayRead(u, &uu);   VecGetArrayRead(v, &vv);   VecGetArrayRead(p, &pp);
  VecGetArrayRead(gradP_x, &gpx); VecGetArrayRead(gradP_y, &gpy);
  VecGetArray(F_0f, &f0f);   VecGetArray(F_f, &ff);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  std::vector<Node*> & node_set = p_mesh->getNodeSet();
  std::map<int, std::vector<Face*> > & face_zone_map = p_mesh->getFaceZoneMap();
  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    if (!is_interior_zone[zone])
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        Face * face = (it->second).at(j);
        long int cell_id1 = face->cell_id1();
        Vec3d face_normal = face->faceNormal();
        double u_face = has_U_BC[zone] ? evaluateFaceU_BC(zone, face, p_mesh) : uu[cell_id1-1];
        double v_face = has_V_BC[zone] ? V_BC[zone] : vv[cell_id1-1];
        double fface = u_face * face_normal.x() + v_face * face_normal.y();
        f0f[face->id()] = fface;

        if (!has_p_BC[zone])
        {
          // Homogeneous Neumann pressure (wall/inlet/slip wall): nothing to correct,
          // matching updateFfaceStar().
          ff[face->id()] = fface;
          continue;
        }

        // pressure_outlet. u_face above is the cell's own solved velocity, which
        // already reflects that cell's momentum balance built with the SMOOTH
        // (cell-centered) pressure gradient (gpx/gpy at cell_id1). Converting to the
        // COMPACT boundary term that M_P/updateFfaceStar() assume (DEB08 3.5.1: normal
        // term only, no cross diffusion) means subtracting the embedded smooth
        // contribution and adding the compact one -- the same "compact minus smooth"
        // logic the interior branch below uses, with the boundary's one-sided compact
        // term standing in for compact_grad_dS. Copying fface straight through, as this
        // branch used to for every non-interior zone, is only correct for the
        // homogeneous-Neumann zones above; see docs/to-do-list-09-22-2026.md item 1.
        const FluentTriCell & cell = cell_set.at(cell_id1-1);
        double distance = p_mesh->node_to_face_distance(cell.centroid(),
                                                         node_set.at(face->node_id1()-1)->point(),
                                                         node_set.at(face->node_id2()-1)->point());
        double compact_dS = face->area() * (p_BC[zone] - pp[cell_id1-1]) / distance;
        double smooth_dS  = gpx[cell_id1-1] * face_normal.x() + gpy[cell_id1-1] * face_normal.y();
        ff[face->id()] = fface - DT * (compact_dS - smooth_dS);
      }
    }
    else
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        Face * face = (it->second).at(j);
        long int face_id = face->id();
        long int cell_id1 = face->cell_id1();
        long int cell_id2 = face->cell_id2();

        double r = face->distance_ratio();
        double u_face = uu[cell_id1-1] * (1.0 - r) + uu[cell_id2-1] * r;
        double v_face = vv[cell_id1-1] * (1.0 - r) + vv[cell_id2-1] * r;

        const FluentTriCell & cell_1 = cell_set.at(cell_id1-1);
        const FluentTriCell & cell_2 = cell_set.at(cell_id2-1);
        double v1 = cell_1.volume();
        double v2 = cell_2.volume();
        const Point & ct1 = cell_1.centroid();
        const Point & ct2 = cell_2.centroid();
        Vec3d ct_to_ct = ct2 - ct1;
        double distance = ct_to_ct.norm();
        Vec3d n1 = ct_to_ct.unitVector();

        Vec3d face_normal = face->faceNormal();
        Vec3d nf = face_normal.unitVector();
        Vec3d n2 = nf - n1;

        double area_f = face->area();
        double alpha_12 = v2 / (v1 + v2);
        double alpha_21 = 1.0 - alpha_12;

        double fface = u_face * face_normal.x() + v_face * face_normal.y();
        f0f[face_id] = fface;

        // Rhie-Chow momentum interpolation. The compact face gradient (DEB-11: the
        // two-point term plus the nonorthogonal correction) is what the pressure
        // equation sees; the smooth part is the cell gradients interpolated to the face.
        // Only their DIFFERENCE belongs here.
        //
        // Subtracting the full compact gradient instead -- as docs/newton_krylov_
        // formulation.tex eq. (5) does, inherited from updateFfaceStar() -- is only
        // correct in the segregated setting, where F_0f is built from the predictor u*.
        // u* has never had the pressure gradient subtracted, so it still carries
        // +dt*interp(grad p) relative to the corrected velocity, which cancels all but
        // this difference. Written monolithically in terms of the velocity that actually
        // solves momentum, that cancellation is absent and the full term is spurious and
        // O(dt) too large.
        double cross_diff = alpha_12 * (gpx[cell_id1-1] * n2.x() + gpy[cell_id1-1] * n2.y())
                          + alpha_21 * (gpx[cell_id2-1] * n2.x() + gpy[cell_id2-1] * n2.y());
        double compact_grad_dS = area_f * ((pp[cell_id2-1] - pp[cell_id1-1]) / distance + cross_diff);

        // face_normal is area-scaled (FluentTwoDMesh.C sets area = |face_normal|), so
        // dotting into it already carries the dS factor.
        double gpx_f = alpha_12 * gpx[cell_id1-1] + alpha_21 * gpx[cell_id2-1];
        double gpy_f = alpha_12 * gpy[cell_id1-1] + alpha_21 * gpy[cell_id2-1];
        double smooth_grad_dS = gpx_f * face_normal.x() + gpy_f * face_normal.y();

        ff[face_id] = fface - DT * (compact_grad_dS - smooth_grad_dS);
      }
    }
  }

  VecRestoreArrayRead(u, &uu);   VecRestoreArrayRead(v, &vv);   VecRestoreArrayRead(p, &pp);
  VecRestoreArrayRead(gradP_x, &gpx); VecRestoreArrayRead(gradP_y, &gpy);
  VecRestoreArray(F_0f, &f0f);   VecRestoreArray(F_f, &ff);
}

// Discrete continuity, sum_f F_f over each cell's faces with the same face orientation
// convention updateMassVeclocities() uses (normal points cell1 -> cell2). Scaled by
// 1/dt to match that function's b_p convention, so the pressure row's magnitude is
// comparable to the momentum rows'. SNES's convergence test takes a single norm over
// the whole residual, so a badly scaled block would silently dominate or vanish.
//
// Written as -(1/dt) sum_f F_f rather than +. The sign of a residual is free -- R = 0
// and -R = 0 are the same equation -- and this choice is what makes M_P directly usable
// as the Schur-complement preconditioner. With +, eliminating velocity gives
// S ~ -L_cmp, the negative of M_P (see NewtonSolverConfig::M_P), so the fieldsplit would
// need a negated copy and GAMG would be handed a negative-definite operator. Flipping
// here flips J_pp and J_pu, hence S, and costs nothing else.
void computeContinuityResidual(FluentTwoDMesh * p_mesh, Vec F_f, Vec p, PetscScalar * div)
{
  const PetscScalar * ff, * pp;
  VecGetArrayRead(F_f, &ff);
  VecGetArrayRead(p, &pp);

  long int n_Cell = (long int)p_mesh->getCellSet().size();
  for (long int i = 0; i < n_Cell; i++)
    div[i] = 0.0;

  std::map<int, std::vector<Face*> > & face_zone_map = p_mesh->getFaceZoneMap();
  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    for (unsigned int j = 0; j < (it->second).size(); j++)
    {
      Face * face = (it->second).at(j);
      double fface = ff[face->id()];
      div[face->cell_id1()-1] -= fface / DT;
      if (is_interior_zone[zone])
        div[face->cell_id2()-1] += fface / DT;
    }
  }

  // Mirrors the MatZeroRows(M_P, row 0) / bb_p[0] = 0.0 pair the segregated path uses:
  // with no pressure B.C. anywhere the pressure is defined only up to a constant, so one
  // continuity row is traded for an anchor. Deliberately left POSITIVE despite the sign
  // flip above, so d(R_p)/dp on this row is +1 and still matches M_P's identity row.
  if (PRESSURE_ANCHOR_NEEDED)
    div[0] = pp[0];

  VecRestoreArrayRead(F_f, &ff);
  VecRestoreArrayRead(p, &pp);
}

PetscErrorCode FormResidual(SNES /*snes*/, Vec X, Vec R, void * ctx)
{
  NewtonCtx * nc = (NewtonCtx *)ctx;
  NewtonSolverConfig * cfg = nc->cfg;
  FluentTwoDMesh * p_mesh = cfg->p_mesh;
  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  long int n_Cell = (long int)cell_set.size();

  // --- unpack X into the per-field scratch vectors -------------------------------
  {
    const PetscScalar * xx;
    PetscScalar * uu, * vv, * pp, * tt;
    VecGetArrayRead(X, &xx);
    VecGetArray(nc->u_w, &uu); VecGetArray(nc->v_w, &vv);
    VecGetArray(nc->p_w, &pp); VecGetArray(nc->T_w, &tt);
    for (long int i = 0; i < n_Cell; i++)
    {
      uu[i] = xx[N_DOF_PER_CELL*i + 0];
      vv[i] = xx[N_DOF_PER_CELL*i + 1];
      pp[i] = xx[N_DOF_PER_CELL*i + 2];
      tt[i] = xx[N_DOF_PER_CELL*i + 3];
    }
    VecRestoreArrayRead(X, &xx);
    VecRestoreArray(nc->u_w, &uu); VecRestoreArray(nc->v_w, &vv);
    VecRestoreArray(nc->p_w, &pp); VecRestoreArray(nc->T_w, &tt);
  }

  // --- face mass flux at the current iterate -------------------------------------
  evaluatePressureGradientValues(p_mesh, nc->p_w, nc->grad_p_x, nc->grad_p_y, cfg->grad_p->data());
  computeFaceFluxes(p_mesh, nc->u_w, nc->v_w, nc->p_w, nc->grad_p_x, nc->grad_p_y, nc->F_0f, nc->F_f);

  // --- momentum: R = M(F_f) * phi - b(F_f) ---------------------------------------
  // The advection operator is rebuilt from F_f, and the RHS from this step's frozen
  // baseline, on every evaluation. Note this includes the pressure-gradient source
  // unconditionally: the segregated Picard loop omits it (DEB-9, "mass velocity") and
  // only adds it in the final re-solve, but the monolithic momentum residual always
  // carries the dV*(dp/dx) term.
  MatCopy(cfg->M_USTAR_FIXED, cfg->M_USTAR, SAME_NONZERO_PATTERN);
  MatCopy(cfg->M_VSTAR_FIXED, cfg->M_VSTAR, SAME_NONZERO_PATTERN);
  VecCopy(nc->b_u_step, nc->b_u_work);
  VecCopy(nc->b_v_step, nc->b_v_work);
  updateAdvectionOperator(p_mesh, nc->F_f, nc->b_u_work, nc->b_v_work, cfg->M_USTAR, cfg->M_VSTAR);
  if (BUOYANCY_ON)
    updateBuoyancySource(p_mesh, nc->T_w, nc->b_v_work);
  updatePressureGradientAsSource(p_mesh, nc->b_u_work, nc->b_v_work, nc->grad_p_x, nc->grad_p_y);

  MatMult(cfg->M_USTAR, nc->u_w, nc->R_u); VecAXPY(nc->R_u, -1.0, nc->b_u_work);
  MatMult(cfg->M_VSTAR, nc->v_w, nc->R_v); VecAXPY(nc->R_v, -1.0, nc->b_v_work);

  // --- continuity -----------------------------------------------------------------
  computeContinuityResidual(p_mesh, nc->F_f, nc->p_w, nc->div.data());

  // --- energy ----------------------------------------------------------------------
  // Same shape as the momentum rows: M_T(F_f) * T - b_T(F_f). Unlike the segregated
  // path, which solves this once per step outside the F_f loop using an already-frozen
  // flux, F_f here is a function of the current iterate, so momentum and energy are
  // converged against each other. Together with the buoyancy source added to b_v_work
  // above -- the B_vT block -- this is what removes the one-time-step temperature lag
  // that makes "-buoyancy true" only a pseudo-transient in the segregated solver.
  MatCopy(cfg->M_T_FIXED, cfg->M_T, SAME_NONZERO_PATTERN);
  VecCopy(nc->b_T_step, nc->b_T_work);
  updateEnergyAdvectionOperator(p_mesh, nc->F_f, nc->b_T_work, cfg->M_T);
  MatMult(cfg->M_T, nc->T_w, nc->R_T); VecAXPY(nc->R_T, -1.0, nc->b_T_work);

  // --- interleave into R ------------------------------------------------------------
  {
    const PetscScalar * ru, * rv, * rt;
    PetscScalar * rr;
    VecGetArrayRead(nc->R_u, &ru); VecGetArrayRead(nc->R_v, &rv); VecGetArrayRead(nc->R_T, &rt);
    VecGetArray(R, &rr);
    for (long int i = 0; i < n_Cell; i++)
    {
      rr[N_DOF_PER_CELL*i + 0] = ru[i];
      rr[N_DOF_PER_CELL*i + 1] = rv[i];
      rr[N_DOF_PER_CELL*i + 2] = nc->div[i];
      rr[N_DOF_PER_CELL*i + 3] = rt[i];
    }
    VecRestoreArrayRead(nc->R_u, &ru); VecRestoreArrayRead(nc->R_v, &rv); VecRestoreArrayRead(nc->R_T, &rt);
    VecRestoreArray(R, &rr);
  }

  return 0;
}

// Two-ring face adjacency, blocked 4x4. This must be a SUPERSET of the true Jacobian
// pattern: too small silently corrupts an FD-coloured Jacobian, too large only costs
// extra colors. Two rings is enough because every dependence path is at most that long
// -- the nonorthogonal cross-diffusion in M_*_FIXED reaches a face neighbour's
// neighbours, and so does F_f through the Green-Gauss grad_p stencil.
void buildJacobianPattern(FluentTwoDMesh * p_mesh, Mat J)
{
  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  long int n_Cell = (long int)cell_set.size();

  std::vector<PetscScalar> zeros;
  std::vector<PetscInt> cols;
  for (long int c = 0; c < n_Cell; c++)
  {
    long int cell_id = cell_set.at(c).id();

    // Boundary faces carry cell_id2 == 0, and FluentTwoDMesh pushes that 0 into the
    // owning cell's neighbour list as a "no neighbour" sentinel, so every id has to be
    // checked before it is used as an index.
    std::set<long int> ring;
    ring.insert(cell_id);
    const std::vector<long int> & nbs = cell_set.at(c).getNeighborCellIDs();
    for (unsigned int i = 0; i < nbs.size(); i++)
    {
      if (nbs[i] <= 0) continue;
      ring.insert(nbs[i]);
      const std::vector<long int> & nbs2 = cell_set.at(nbs[i]-1).getNeighborCellIDs();
      for (unsigned int j = 0; j < nbs2.size(); j++)
        if (nbs2[j] > 0)
          ring.insert(nbs2[j]);
    }

    cols.clear();
    for (std::set<long int>::const_iterator it = ring.begin(); it != ring.end(); ++it)
      for (int d = 0; d < N_DOF_PER_CELL; d++)
        cols.push_back(N_DOF_PER_CELL * (*it - 1) + d);

    zeros.assign(cols.size(), 0.0);
    for (int d = 0; d < N_DOF_PER_CELL; d++)
    {
      PetscInt row = N_DOF_PER_CELL * (cell_id - 1) + d;
      MatSetValues(J, 1, &row, (PetscInt)cols.size(), cols.data(), zeros.data(), INSERT_VALUES);
    }
  }

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  // Keep the explicit zeros: they ARE the pattern the colouring is derived from.
  MatSetOption(J, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
}

} // anonymous namespace

void runNewtonSolver(NewtonSolverConfig & cfg)
{
  FluentTwoDMesh * p_mesh = cfg.p_mesh;
  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  long int n_Cell = (long int)cell_set.size();
  PetscInt n_dof = (PetscInt)(N_DOF_PER_CELL * n_Cell);

  NewtonCtx nc;
  nc.cfg = &cfg;

  VecDuplicate(cfg.u, &nc.u_w);       VecDuplicate(cfg.u, &nc.v_w);
  VecDuplicate(cfg.u, &nc.p_w);       VecDuplicate(cfg.u, &nc.T_w);
  VecDuplicate(cfg.u, &nc.grad_p_x);  VecDuplicate(cfg.u, &nc.grad_p_y);
  VecDuplicate(cfg.u, &nc.b_u_work);  VecDuplicate(cfg.u, &nc.b_v_work);
  VecDuplicate(cfg.u, &nc.b_T_work);
  VecDuplicate(cfg.u, &nc.b_u_step);  VecDuplicate(cfg.u, &nc.b_v_step);
  VecDuplicate(cfg.u, &nc.b_T_step);
  VecDuplicate(cfg.u, &nc.b_u_bc);    VecDuplicate(cfg.u, &nc.b_v_bc);
  VecDuplicate(cfg.u, &nc.b_T_bc);
  VecDuplicate(cfg.u, &nc.u_old);     VecDuplicate(cfg.u, &nc.v_old);
  VecDuplicate(cfg.u, &nc.T_old);
  VecDuplicate(cfg.u, &nc.R_u);       VecDuplicate(cfg.u, &nc.R_v);
  VecDuplicate(cfg.u, &nc.R_T);
  VecDuplicate(cfg.F_face, &nc.F_0f); VecDuplicate(cfg.F_face, &nc.F_f);
  nc.div.resize(n_Cell);

  // Recover the time-independent part of the RHS baselines. On entry cxx_main.C's
  // b_USTAR holds (B.C. diffusion contributions) + u^0 * rho*V/dt, and cfg.u still holds
  // u^0, so subtracting the latter leaves the B.C. part alone. Doing this once here and
  // rebuilding b = b_bc + phi^n * rho*V/dt per step makes "phi^n is frozen across the
  // Newton iterations" structurally true, rather than a property maintained by the
  // in-place roll-forward the segregated loop relies on.
  {
    const PetscScalar * uu, * vv, * tt, * bu, * bv, * bt;
    PetscScalar * bu_bc, * bv_bc, * bt_bc;
    VecGetArrayRead(cfg.u, &uu); VecGetArrayRead(cfg.v, &vv); VecGetArrayRead(cfg.T, &tt);
    VecGetArrayRead(cfg.b_USTAR, &bu); VecGetArrayRead(cfg.b_VSTAR, &bv); VecGetArrayRead(cfg.b_T, &bt);
    VecGetArray(nc.b_u_bc, &bu_bc); VecGetArray(nc.b_v_bc, &bv_bc); VecGetArray(nc.b_T_bc, &bt_bc);
    const double RHO = 1.0;
    for (long int i = 0; i < n_Cell; i++)
    {
      double vol_dt = RHO * cell_set.at(i).volume() / DT;
      bu_bc[i] = bu[i] - uu[i] * vol_dt;
      bv_bc[i] = bv[i] - vv[i] * vol_dt;
      bt_bc[i] = bt[i] - tt[i] * vol_dt;
    }
    VecRestoreArrayRead(cfg.u, &uu); VecRestoreArrayRead(cfg.v, &vv); VecRestoreArrayRead(cfg.T, &tt);
    VecRestoreArrayRead(cfg.b_USTAR, &bu); VecRestoreArrayRead(cfg.b_VSTAR, &bv); VecRestoreArrayRead(cfg.b_T, &bt);
    VecRestoreArray(nc.b_u_bc, &bu_bc); VecRestoreArray(nc.b_v_bc, &bv_bc); VecRestoreArray(nc.b_T_bc, &bt_bc);
  }

  // --- state vector, residual, Jacobian -------------------------------------------
  Vec X, R;
  VecCreate(PETSC_COMM_SELF, &X);
  VecSetSizes(X, PETSC_DECIDE, n_dof);
  VecSetFromOptions(X);
  VecDuplicate(X, &R);

  Mat J;
  // Two rings of a triangle is at most 1 + 3 + 3*2 = 10 cells, so 40 columns per row;
  // 64 leaves headroom for irregular connectivity without a reallocation.
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_dof, n_dof, 64, NULL, &J);
  // Lets PCFIELDSPLIT recover the four fields from the interlaced layout, so a
  // Schur-complement split can be driven entirely from the command line.
  MatSetBlockSize(J, N_DOF_PER_CELL);
  MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  buildJacobianPattern(p_mesh, J);

  SNES snes;
  SNESCreate(PETSC_COMM_SELF, &snes);
  SNESSetFunction(snes, R, FormResidual, &nc);

  // Finite-difference the Jacobian, coloured from J's nonzero pattern. This costs one
  // residual evaluation per colour, which is slow, but it puts all the correctness risk
  // in FormResidual() -- if the answer is wrong, the residual is wrong. Swapping in
  // -snes_mf_operator with a physics-based preconditioner later reuses the same residual.
  // Passing a NULL colouring context lets SNESComputeJacobianDefaultColor() build the
  // MatFDColoring itself from J's pattern on first use, which is why buildJacobianPattern()
  // above has to get that pattern right.
  SNESSetJacobian(snes, J, J, SNESComputeJacobianDefaultColor, NULL);

  // Default to a direct solve. The Rhie-Chow pressure operator is the difference between
  // the compact and the reconstructed Laplacian -- nearly null on smooth modes, by
  // design, since its job is to see checkerboard modes and little else -- so the coupled
  // system is a genuinely stiff saddle point and PETSc's default GMRES+ILU(0) diverges
  // on it at the first Newton step. This solver is serial throughout (PETSC_COMM_SELF)
  // and the largest shipped mesh is ~21k cells, so LU is affordable and robust. Set
  // before SNESSetFromOptions() so -ksp_type/-pc_type on the command line still win;
  // swapping in a physics-based preconditioner later is the intended upgrade path.
  {
    KSP ksp; PC pc;
    SNESGetKSP(snes, &ksp);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
  }
  SNESSetFromOptions(snes);

  // With -newton_schur_mp, hand M_P to a Schur fieldsplit as the preconditioner for the
  // pressure block. PETSc's built-in "selfp" instead builds D diag(A)^-1 G -- the WIDE
  // -stencil Laplacian scaled by dV/a_P, the rough-mode response -- where the derivation
  // in newton_solver.h says the answer is the COMPACT one at full strength. Opt-in so
  // that -pc_fieldsplit_schur_precondition still means what it says when not passed, and
  // after SNESSetFromOptions(), which is what creates the PCFIELDSPLIT.
  {
    PetscBool use_mp = PETSC_FALSE, is_fs = PETSC_FALSE;
    PetscOptionsGetBool(NULL, NULL, "-newton_schur_mp", &use_mp, NULL);
    KSP ksp; PC pc;
    SNESGetKSP(snes, &ksp);
    KSPGetPC(ksp, &pc);
    PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_fs);
    if (use_mp && is_fs && cfg.M_P)
      PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, cfg.M_P);
  }

  // --- time loop -------------------------------------------------------------------
  // Output mirrors the segregated loop exactly (same records, same interval rule, same
  // final-step stream function) so the two modes' Exodus files are directly comparable.
  int output_index = 0;
  const double RHO = 1.0;
  for (int time_step = 0; time_step < cfg.n_timestep; time_step++)
  {
    std::cout << "Time Step " << time_step << " ... ";

    VecCopy(cfg.u, nc.u_old);
    VecCopy(cfg.v, nc.v_old);
    VecCopy(cfg.T, nc.T_old);

    // b_*_step = b_*_bc + phi^n * rho*V/dt
    {
      const PetscScalar * uu, * vv, * tt, * bu_bc, * bv_bc, * bt_bc;
      PetscScalar * bu, * bv, * bt;
      VecGetArrayRead(nc.u_old, &uu); VecGetArrayRead(nc.v_old, &vv); VecGetArrayRead(nc.T_old, &tt);
      VecGetArrayRead(nc.b_u_bc, &bu_bc); VecGetArrayRead(nc.b_v_bc, &bv_bc); VecGetArrayRead(nc.b_T_bc, &bt_bc);
      VecGetArray(nc.b_u_step, &bu); VecGetArray(nc.b_v_step, &bv); VecGetArray(nc.b_T_step, &bt);
      for (long int i = 0; i < n_Cell; i++)
      {
        double vol_dt = RHO * cell_set.at(i).volume() / DT;
        bu[i] = bu_bc[i] + uu[i] * vol_dt;
        bv[i] = bv_bc[i] + vv[i] * vol_dt;
        bt[i] = bt_bc[i] + tt[i] * vol_dt;
      }
      VecRestoreArrayRead(nc.u_old, &uu); VecRestoreArrayRead(nc.v_old, &vv); VecRestoreArrayRead(nc.T_old, &tt);
      VecRestoreArrayRead(nc.b_u_bc, &bu_bc); VecRestoreArrayRead(nc.b_v_bc, &bv_bc); VecRestoreArrayRead(nc.b_T_bc, &bt_bc);
      VecRestoreArray(nc.b_u_step, &bu); VecRestoreArray(nc.b_v_step, &bv); VecRestoreArray(nc.b_T_step, &bt);
    }

    // Initial guess: the previous time level.
    {
      const PetscScalar * uu, * vv, * pp, * tt;
      PetscScalar * xx;
      VecGetArrayRead(cfg.u, &uu); VecGetArrayRead(cfg.v, &vv);
      VecGetArrayRead(cfg.p, &pp); VecGetArrayRead(cfg.T, &tt);
      VecGetArray(X, &xx);
      for (long int i = 0; i < n_Cell; i++)
      {
        xx[N_DOF_PER_CELL*i + 0] = uu[i];
        xx[N_DOF_PER_CELL*i + 1] = vv[i];
        xx[N_DOF_PER_CELL*i + 2] = pp[i];
        xx[N_DOF_PER_CELL*i + 3] = tt[i];
      }
      VecRestoreArrayRead(cfg.u, &uu); VecRestoreArrayRead(cfg.v, &vv);
      VecRestoreArrayRead(cfg.p, &pp); VecRestoreArrayRead(cfg.T, &tt);
      VecRestoreArray(X, &xx);
    }

    SNESSolve(snes, NULL, X);

    SNESConvergedReason reason;
    SNESGetConvergedReason(snes, &reason);
    if (reason < 0)
    {
      PetscInt its; SNESGetIterationNumber(snes, &its);
      std::cout << std::endl << "ERROR: SNES diverged at time step " << time_step
                << " (reason " << (int)reason << ", " << (int)its << " iterations)."
                << std::endl;
      std::cout << "       Re-run with -snes_monitor -snes_converged_reason for detail."
                << std::endl;
      exit(1);
    }

    // Unpack the converged state.
    {
      const PetscScalar * xx;
      PetscScalar * uu, * vv, * pp, * tt;
      VecGetArrayRead(X, &xx);
      VecGetArray(cfg.u, &uu); VecGetArray(cfg.v, &vv);
      VecGetArray(cfg.p, &pp); VecGetArray(cfg.T, &tt);
      for (long int i = 0; i < n_Cell; i++)
      {
        uu[i] = xx[N_DOF_PER_CELL*i + 0];
        vv[i] = xx[N_DOF_PER_CELL*i + 1];
        pp[i] = xx[N_DOF_PER_CELL*i + 2];
        tt[i] = xx[N_DOF_PER_CELL*i + 3];
      }
      VecRestoreArrayRead(X, &xx);
      VecRestoreArray(cfg.u, &uu); VecRestoreArray(cfg.v, &vv);
      VecRestoreArray(cfg.p, &pp); VecRestoreArray(cfg.T, &tt);
    }

    // Face flux at the converged state, for the stream function on the final step.
    evaluatePressureGradientValues(p_mesh, cfg.p, nc.grad_p_x, nc.grad_p_y, cfg.grad_p->data());
    computeFaceFluxes(p_mesh, cfg.u, cfg.v, cfg.p, nc.grad_p_x, nc.grad_p_y, nc.F_0f, cfg.F_face);

    // Per-step change, same definition and same CSV columns as the segregated loop.
    PetscScalar error_norm_u, error_norm_v, error_norm_T;
    VecAXPY(nc.u_old, -1.0, cfg.u); VecNorm(nc.u_old, NORM_2, &error_norm_u);
    VecAXPY(nc.v_old, -1.0, cfg.v); VecNorm(nc.v_old, NORM_2, &error_norm_v);
    VecAXPY(nc.T_old, -1.0, cfg.T); VecNorm(nc.T_old, NORM_2, &error_norm_T);
    *(cfg.convergence_file) << time_step << "," << error_norm_u << "," << error_norm_v
                            << "," << error_norm_T << std::endl;

    bool computeStreamFunction = (time_step == cfg.n_timestep - 1);
    if ((time_step + 1) % cfg.output_interval == 0 || computeStreamFunction)
    {
      writeExodusStep(cfg.exoCtx, p_mesh, output_index, (time_step + 1) * DT,
                      cfg.u, cfg.v, cfg.p, cfg.T, cfg.F_face, computeStreamFunction);
      output_index++;
    }

    std::cout << "  [ok]" << std::endl;
  }

  SNESDestroy(&snes);
  MatDestroy(&J);
  VecDestroy(&X); VecDestroy(&R);
  VecDestroy(&nc.u_w); VecDestroy(&nc.v_w); VecDestroy(&nc.p_w); VecDestroy(&nc.T_w);
  VecDestroy(&nc.grad_p_x); VecDestroy(&nc.grad_p_y);
  VecDestroy(&nc.F_0f); VecDestroy(&nc.F_f);
  VecDestroy(&nc.b_u_work); VecDestroy(&nc.b_v_work); VecDestroy(&nc.b_T_work);
  VecDestroy(&nc.b_u_step); VecDestroy(&nc.b_v_step); VecDestroy(&nc.b_T_step);
  VecDestroy(&nc.b_u_bc); VecDestroy(&nc.b_v_bc); VecDestroy(&nc.b_T_bc);
  VecDestroy(&nc.u_old); VecDestroy(&nc.v_old); VecDestroy(&nc.T_old);
  VecDestroy(&nc.R_u); VecDestroy(&nc.R_v); VecDestroy(&nc.R_T);
}
