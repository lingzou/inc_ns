#include <fstream>
#include <iostream>
#include <vector>
#include <petscsnes.h>

#include "inc_problems.h"
#include "ExodusOutput.h"
#include "newton_solver.h"

/*
 * References:
 * 1. Amaresh Dalal, V. Eswaran, and G. Biswas, "A Finite-Volume Method for Navier-Stokes Equations on Unstructured Meshes,"
 *    Numerical Heat Transfer, Part B, 54: 238-259, 2008.
 */

struct PETSC_APPCTX
{
  FluentTwoDMesh * p_mesh;

  void FreeWorkSpace();
};

void PETSC_APPCTX::FreeWorkSpace()
{
  if (p_mesh != NULL)
    delete p_mesh;
}


int main(int argc, char **argv)
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  char      inputFile[PETSC_MAX_PATH_LEN];
  PetscBool input_specified = PETSC_FALSE;
  PetscOptionsGetString(NULL, NULL, "-i", inputFile, sizeof(inputFile), &input_specified);

  if (!input_specified)
  {
    std::cout << "ERROR: Input file is not specified." << std::endl;
    std::cout << "Please specify the input file using '-i <input_file>', for example, '-i cylinder.i'." << std::endl;
    exit(1);
  }

  // PETSC_TRUE ("require") makes PETSc itself raise a clear error if inputFile
  // doesn't exist, so no separate existence check is needed here.
  PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, inputFile, PETSC_TRUE);

  char      problem_name[PETSC_MAX_PATH_LEN];
  char      mesh_file[PETSC_MAX_PATH_LEN];
  PetscBool problem_specified = PETSC_FALSE;
  PetscBool mesh_specified = PETSC_FALSE;
  PetscBool dt_specified = PETSC_FALSE;
  PetscBool visc_specified = PETSC_FALSE;
  PetscBool n_timestep_specified = PETSC_FALSE;
  PetscInt  n_timestep_opt = 0;
  PetscInt  output_interval = 1;
  PetscBool buoyancy_flg = PETSC_FALSE;
  PetscBool newton_flg = PETSC_FALSE;
  double    T_init = 0.0;

  PetscOptionsGetString(NULL, NULL, "-problem", problem_name, sizeof(problem_name), &problem_specified);
  PetscOptionsGetString(NULL, NULL, "-mesh", mesh_file, sizeof(mesh_file), &mesh_specified);
  PetscOptionsGetReal(NULL, NULL, "-dt", &DT, &dt_specified);
  PetscOptionsGetReal(NULL, NULL, "-visc", &VISC, &visc_specified);
  // -diff_t is optional (default 0.0): problems with no T Dirichlet B.C. anywhere keep
  // T uniformly at -T_init forever regardless of its value, so existing input files
  // that predate the energy equation don't need to specify it.
  PetscOptionsGetReal(NULL, NULL, "-diff_t", &DIFF_T, NULL);
  PetscOptionsGetInt(NULL, NULL, "-n_timestep", &n_timestep_opt, &n_timestep_specified);
  // -output_interval is optional (default 1 = write every time step, the historical
  // behavior). Long pseudo-transient runs to steady state only need a coarse history:
  // cavity_buoyancy_1e5's 4000 steps write 1.1 GB at interval 1.
  PetscOptionsGetInt(NULL, NULL, "-output_interval", &output_interval, NULL);
  PetscOptionsGetBool(NULL, NULL, "-buoyancy", &buoyancy_flg, NULL);
  BUOYANCY_ON = (bool)buoyancy_flg;
  PetscOptionsGetReal(NULL, NULL, "-T_init", &T_init, NULL);
  // -newton is optional (default false = the segregated Picard/SIMPLE path below).
  // Read here, after the PetscOptionsInsertFile() above, so it can be set either on the
  // command line or from the -i input file like every other option.
  PetscOptionsGetBool(NULL, NULL, "-newton", &newton_flg, NULL);
  NEWTON_ON = (bool)newton_flg;

  if (!problem_specified || !mesh_specified || !dt_specified || !visc_specified || !n_timestep_specified)
  {
    std::cout << "ERROR: '" << inputFile << "' must specify all of: -problem, -mesh, -dt, -visc, -n_timestep." << std::endl;
    exit(1);
  }

  if (output_interval < 1)
  {
    std::cout << "ERROR: '-output_interval' must be >= 1 (got " << output_interval << ")." << std::endl;
    exit(1);
  }

  PETSC_APPCTX app;

  app.p_mesh = new FluentTwoDMesh();
  app.p_mesh->createMeshFromFile(mesh_file, false, true);

  // Per-zone B.C.s (-bc_<zone>_type ...) are read once the mesh's zone IDs are known,
  // so every boundary zone the mesh actually has can be checked against the input file.
  setupBoundaryConditions(app.p_mesh);

  // Derive the output basename from the input file, e.g. "cylinder.i" -> "cylinder", so
  // running with different -i files never overwrites each other's output.
  std::string input_path(inputFile);
  size_t slash = input_path.find_last_of("/\\");
  std::string input_base = (slash == std::string::npos) ? input_path : input_path.substr(slash + 1);
  size_t dot = input_base.find_last_of('.');
  if (dot != std::string::npos)
    input_base = input_base.substr(0, dot);
  std::string output_prefix = "output/" + input_base + "_out";
  // The Newton path writes to its own "<base>_out_newton.*" so that running the same
  // input file both ways leaves two results side by side to compare -- the refinement
  // -rate check in runNewtonTests.sh depends on this -- and so that -newton can never
  // silently overwrite an expected/ gold produced by the segregated solver.
  if (NEWTON_ON)
    output_prefix += "_newton";

  ExodusOutputContext * exoCtx = createExodusOutput(app.p_mesh, (output_prefix + ".e").c_str());

  // Per-time-step convergence history (u/v L2 error), flushed after every row so an
  // external process (e.g. gui.py) can tail it live while this run is still going.
  std::ofstream convergence_file(output_prefix + "_convergence.csv");
  convergence_file << "time_step,u_L2_error,v_L2_error,T_L2_error" << std::endl;

  std::cout << "# of faces: " << app.p_mesh->n_Faces() << std::endl;

  std::vector<Node*> & node_set = app.p_mesh->getNodeSet();
  std::map<int, std::vector<Face*> > & face_zone_map = app.p_mesh->getFaceZoneMap();
  std::vector<FluentTriCell> & cell_set = app.p_mesh->getCellSet();

  std::cout << "Number of nodes = " << node_set.size() << std::endl;
  std::cout << "Number of cells = " << cell_set.size() << std::endl;

  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    long int zone_id = it->first;
    std::vector<Face*> & faces = it->second;

    std::cout << "Face Zone " << zone_id << std::endl;
    std::cout << "   number of faces = " << faces.size() << std::endl;
  }

  long int n_Cell = app.p_mesh->n_Cells();
  long int n_Face = app.p_mesh->n_Faces();

  // --->
  std::cout << "GRAD START" << std::endl;
  // Heap-allocated (not stack arrays): n_Cell scales with mesh size, and
  // backward_facing_step's ~21k cells * 4 GRAD arrays already exceeds a typical 8MB
  // thread stack.
  std::vector<GRAD> grad_u_star(n_Cell);
  std::vector<GRAD> grad_v_star(n_Cell);
  std::vector<GRAD> grad_p(n_Cell);
  std::vector<GRAD> grad_T(n_Cell);
  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    std::cout << "zone = " << zone << std::endl;
    if (!is_interior_zone[zone])
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        Face * face = (it->second)[j];
        long int cell_id1 = face->cell_id1();
        Vec3d face_normal = face->faceNormal();
        double v1 = cell_set.at(cell_id1-1).volume();

        // \grad phi = 1/v * \Sigma phi_face \vec{A_face}; equation on page 244 of Ref. [1], following equation (10).
        // If has boundary value, it goes to the BC contribution; otherwise,
        // we assuming d_phi/dn = 0, such that phi_BC = phi_cell, and this boundary contributes as a cell value
        if (has_U_BC[zone])
          grad_u_star[cell_id1-1].addBCContribution(U_BC[zone] * face_normal.x() / v1, U_BC[zone] * face_normal.y() / v1);
        else // u face = u cell
          grad_u_star[cell_id1-1].addCoef(0, cell_id1, face_normal.x() / v1, face_normal.y() / v1);

        if (has_V_BC[zone])
          grad_v_star[cell_id1-1].addBCContribution(V_BC[zone] * face_normal.x() / v1, V_BC[zone] * face_normal.y() / v1);
        else
          grad_v_star[cell_id1-1].addCoef(0, cell_id1, face_normal.x() / v1, face_normal.y() / v1);

        if (has_p_BC[zone])
          grad_p[cell_id1-1].addBCContribution(p_BC[zone] * face_normal.x() / v1, p_BC[zone] * face_normal.y() / v1);
        else
          grad_p[cell_id1-1].addCoef(0, cell_id1, face_normal.x() / v1, face_normal.y() / v1);

        if (has_T_BC[zone])
          grad_T[cell_id1-1].addBCContribution(T_BC[zone] * face_normal.x() / v1, T_BC[zone] * face_normal.y() / v1);
        else // adiabatic: T face = T cell
          grad_T[cell_id1-1].addCoef(0, cell_id1, face_normal.x() / v1, face_normal.y() / v1);
      }
    }
    else
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        Face * face = (it->second)[j];
        long int cell_id1 = face->cell_id1();
        long int cell_id2 = face->cell_id2();
        double v1 = cell_set.at(cell_id1-1).volume();
        double v2 = cell_set.at(cell_id2-1).volume();
        Vec3d face_normal = face->faceNormal();  // face norm is always pointing from cell 1 to cell 2, this is how mesh was prepared

        // phi_face = (v2 * phi_1 + v1 * phi_2) / (v1 + v2); equation (5) of ref. [1]

        double alpha_12 = v2 / (v1 + v2);
        double alpha_21 = 1.0 - alpha_12;
        // Cell 1
        grad_u_star[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);
        grad_v_star[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);
             grad_p[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);
             grad_T[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);

        grad_u_star[cell_id1-1].size++;
        grad_v_star[cell_id1-1].size++;
             grad_p[cell_id1-1].size++;
             grad_T[cell_id1-1].size++;

        grad_u_star[cell_id1-1].addCoef(grad_u_star[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
        grad_v_star[cell_id1-1].addCoef(grad_v_star[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
             grad_p[cell_id1-1].addCoef(     grad_p[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
             grad_T[cell_id1-1].addCoef(     grad_T[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
        // Cell 2
        grad_u_star[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);
        grad_v_star[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);
             grad_p[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);
             grad_T[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);

        grad_u_star[cell_id2-1].size++;
        grad_v_star[cell_id2-1].size++;
             grad_p[cell_id2-1].size++;
             grad_T[cell_id2-1].size++;

        grad_u_star[cell_id2-1].addCoef(grad_u_star[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
        grad_v_star[cell_id2-1].addCoef(grad_v_star[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
             grad_p[cell_id2-1].addCoef(     grad_p[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
             grad_T[cell_id2-1].addCoef(     grad_T[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
      }
    }
  }
  std::cout << "GRAD END" << std::endl;
  // <---

  // KSP stuffs
  std::cout << "KSP START" << std::endl;
  Mat       M_USTAR_FIXED, M_VSTAR_FIXED;
  Mat       M_USTAR, M_VSTAR;
  Mat       M_P;
  Mat       M_T_FIXED, M_T;
  Vec       u, u_STAR, b_USTAR, b_USTAR_FINAL;
  Vec       v, v_STAR, b_VSTAR, b_VSTAR_FINAL;
  Vec       p, b_p, p_old_it, b_p_FINAL;
  Vec       T, T_NEW, b_T, b_T_FINAL;
  Vec       F_face_star, F_0f_star, F_face_star_old_it;
  KSP       ksp_USTAR, ksp_VSTAR, ksp_P, ksp_T;

  Vec       grad_p_x, grad_p_y;
  Vec       p_src_x, p_src_y;

  std::cout << "KSP START SETUP" << std::endl;

  VecCreate(PETSC_COMM_SELF, &u_STAR);
  VecSetSizes(u_STAR, PETSC_DECIDE, n_Cell);
  VecSetFromOptions(u_STAR);
  VecDuplicate(u_STAR, &b_USTAR);
  VecDuplicate(u_STAR, &v_STAR);
  VecDuplicate(u_STAR, &b_VSTAR);
  VecDuplicate(u_STAR, &b_USTAR_FINAL);
  VecDuplicate(u_STAR, &b_VSTAR_FINAL);
  VecDuplicate(u_STAR, &u);
  VecDuplicate(u_STAR, &v);
  VecDuplicate(u_STAR, &p);
  VecDuplicate(u_STAR, &b_p);
  VecDuplicate(u_STAR, &b_p_FINAL);
  VecDuplicate(u_STAR, &p_old_it);

  VecDuplicate(u_STAR, &T);
  VecDuplicate(u_STAR, &T_NEW);
  VecDuplicate(u_STAR, &b_T);
  VecDuplicate(u_STAR, &b_T_FINAL);

  VecDuplicate(u_STAR, &grad_p_x);
  VecDuplicate(u_STAR, &grad_p_y);
  VecDuplicate(u_STAR, &p_src_x);
  VecDuplicate(u_STAR, &p_src_y);

  VecSet(u_STAR, 1.0);
  VecSet(v_STAR, 0.0);
  VecSet(b_USTAR, 0.0);
  VecSet(b_VSTAR, 0.0);
  VecSet(b_USTAR_FINAL, 0.0);
  VecSet(b_VSTAR_FINAL, 0.0);
  VecSet(u, 1.0);
  VecSet(v, 0.0);
  VecSet(p, 0.0);
  VecSet(b_p, 0.0);
  VecSet(b_p_FINAL, 0.0);
  VecSet(p_old_it, 0.0);

  VecSet(T, T_init);
  VecSet(T_NEW, T_init);
  VecSet(b_T, 0.0);
  VecSet(b_T_FINAL, 0.0);

  VecSet(grad_p_x, 0.0);
  VecSet(grad_p_y, 0.0);
  VecSet(p_src_x, 0.0);
  VecSet(p_src_y, 0.0);

  VecCreate(PETSC_COMM_SELF, &F_face_star);
  VecSetSizes(F_face_star, PETSC_DECIDE, n_Face);
  VecSetFromOptions(F_face_star);
  VecDuplicate(F_face_star, &F_0f_star);
  VecDuplicate(F_face_star, &F_face_star_old_it);
  VecSet(F_face_star, 0.0);
  VecSet(F_0f_star, 0.0);
  VecSet(F_face_star_old_it, 0.0);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_USTAR_FIXED);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_USTAR_FIXED, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_VSTAR_FIXED);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_VSTAR_FIXED, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_P);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_P, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_T_FIXED);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_T_FIXED, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // M_USTAR/M_VSTAR are created once (via MatDuplicate from the assembled *_FIXED
  // matrices, once those are finalized below) and reused for the whole run -- each
  // Picard iteration resets their values with MatCopy() rather than MatDuplicate()'ing
  // a brand new matrix every time, which used to leak one Mat per iteration (up to
  // ~100/timestep, never MatDestroy()'d) and paid full allocation cost on top.
  KSPCreate(PETSC_COMM_WORLD, &ksp_USTAR);
  KSPCreate(PETSC_COMM_WORLD, &ksp_VSTAR);
  KSPCreate(PETSC_COMM_WORLD, &ksp_P);
  KSPCreate(PETSC_COMM_WORLD, &ksp_T);
  KSPSetType(ksp_P, KSPCG);
  // M_P discretizes a Poisson-like pressure equation, so PCGAMG is probably the best preconditioner
  PC pc_P;
  KSPGetPC(ksp_P, &pc_P);
  PCSetType(pc_P, PCGAMG);
  // Distinct options prefixes let each solve's -ksp_type/-pc_type/etc. be tuned
  // independently from the command line or input file (e.g. "-p_pc_type gamg"),
  // instead of one shared, unprefixed options namespace covering all three KSPs.
  KSPSetOptionsPrefix(ksp_USTAR, "ustar_");
  KSPSetOptionsPrefix(ksp_VSTAR, "vstar_");
  KSPSetOptionsPrefix(ksp_P, "p_");
  KSPSetOptionsPrefix(ksp_T, "t_");

  std::cout << "KSP END SETUP" << std::endl;

  PetscScalar * bb_ustar;
  PetscScalar * bb_vstar;
  PetscScalar * bb_p;
  PetscScalar * bb_t;
  PetscScalar * uu; //, * uu_star;
  PetscScalar * vv; //, * vv_star;
  PetscScalar * tt;
  VecGetArray(b_USTAR, &bb_ustar);
  VecGetArray(b_VSTAR, &bb_vstar);
  VecGetArray(b_p, &bb_p);
  VecGetArray(b_T, &bb_t);
  VecGetArray(u, &uu);
  VecGetArray(v, &vv);
  VecGetArray(T, &tt);
  //VecGetArray(u_STAR, &uu_star);
  //VecGetArray(v_STAR, &vv_star);

  // Diffusion terms
  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    //double VISC = 1.0;
    int zone = it->first;
    if (!is_interior_zone[zone])
    {
      std::cout << "zone = " << zone << std::endl;
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        Face * face = (it->second).at(j);
        long int cell_id1 = face->cell_id1();
        long int cell_id2 = face->cell_id2();

        long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
        const FluentTriCell & cell = cell_set.at(cell_id-1);

        long int node_id1 = face->node_id1();
        long int node_id2 = face->node_id2();

        Point pt0 = cell.centroid();
        Point pt1 = node_set.at(node_id1-1)->point();
        Point pt2 = node_set.at(node_id2-1)->point();
        double distance = app.p_mesh->node_to_face_distance(pt0, pt1, pt2);

        PetscInt row = cell_id - 1;
        PetscInt col[1]; col[0] = row;
        double val[1];
        val[0] = VISC * face->area() / distance;

        if (has_U_BC[zone])
        {
          MatSetValues(M_USTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
          bb_ustar[row] += VISC * face->area() / distance * U_BC[zone];
        }
        if (has_V_BC[zone])
        {
          MatSetValues(M_VSTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
          bb_vstar[row] += VISC * face->area() / distance * V_BC[zone];
        }
        if (has_p_BC[zone])
        {
          val[0] = -face->area() / distance;
          MatSetValues(M_P, 1, &row, 1, col, val, ADD_VALUES);
          bb_p[row]     -= face->area() / distance * p_BC[zone];
        }
        if (has_T_BC[zone])
        {
          val[0] = DIFF_T * face->area() / distance;
          MatSetValues(M_T_FIXED, 1, &row, 1, col, val, ADD_VALUES);
          bb_t[row] += DIFF_T * face->area() / distance * T_BC[zone];
        }
      }
    }
    else
    {
      std::cout << "zone = " << zone << std::endl;
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
          Face * face = (it->second).at(j);
          long int cell_id1 = face->cell_id1();
          long int cell_id2 = face->cell_id2();

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
          // cell 1
          PetscInt r1 = cell_id1 - 1; PetscInt r2 = cell_id2 - 1;
          PetscInt col[1];
          double val[1], nval[1];

          col[0] = r1; val[0] = VISC * area_f / distance; nval[0] = -val[0];
          MatSetValues(M_USTAR_FIXED, 1, &r1, 1, col, val, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, 1, col, nval, ADD_VALUES);  // cell 2 off-diag
          // copy & paste for VSTAR and (-1.0)MP
          MatSetValues(M_VSTAR_FIXED, 1, &r1, 1, col, val, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, 1, col, nval, ADD_VALUES);  // cell 2 off-diag
          val[0] = area_f / distance; nval[0] = -val[0];
          MatSetValues(M_P, 1, &r1, 1, col, nval, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_P, 1, &r2, 1, col, val, ADD_VALUES);  // cell 2 off-diag
          val[0] = DIFF_T * area_f / distance; nval[0] = -val[0];
          MatSetValues(M_T_FIXED, 1, &r1, 1, col, val, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_T_FIXED, 1, &r2, 1, col, nval, ADD_VALUES);  // cell 2 off-diag

          col[0] = r2; val[0] = VISC * area_f / distance; nval[0] = -val[0];
          MatSetValues(M_USTAR_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag
          // copy & paste for VSTAR and (-1.0)MP
          MatSetValues(M_VSTAR_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag
          val[0] = area_f / distance; nval[0] = -val[0];
          MatSetValues(M_P, 1, &r1, 1, col, val, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_P, 1, &r2, 1, col, nval, ADD_VALUES);   // cell 2 diag
          val[0] = DIFF_T * area_f / distance; nval[0] = -val[0];
          MatSetValues(M_T_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_T_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag

          int size1 = grad_u_star[cell_id1-1].size;
          if (size1 > 3) {std::cerr<<"size1 > 3\n"; exit(1);}
          PetscInt col1[size1+1]; double val1[size1+1], nval1[size1+1];
          for (int i = 0; i < size1+1; i++)
          {
            col1[i] = grad_u_star[cell_id1-1].cell_id[i]-1;
            val1[i] = -grad_u_star[cell_id1-1].coef_x[i] * n2.x() - grad_u_star[cell_id1-1].coef_y[i] * n2.y();
            val1[i] *= (VISC * area_f * alpha_12);
            nval1[i] = -val1[i];
          }
          MatSetValues(M_USTAR_FIXED, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          double bc_contribution = VISC * area_f * alpha_12 * (grad_u_star[cell_id1-1].bc_x * n2.x() + grad_u_star[cell_id1-1].bc_y * n2.y());
          bb_ustar[r1] += bc_contribution;
          bb_ustar[r2] -= bc_contribution;
          // copy & paste for VSTAR
          if (grad_v_star[cell_id1-1].size != size1) {std::cerr<<"size != size1\n"; exit(1);}
          for (int i = 0; i < size1+1; i++)
          {
            col1[i] = grad_v_star[cell_id1-1].cell_id[i]-1;
            val1[i] = -grad_v_star[cell_id1-1].coef_x[i] * n2.x() - grad_v_star[cell_id1-1].coef_y[i] * n2.y();
            val1[i] *= (VISC * area_f * alpha_12);
            nval1[i] = -val1[i];
          }
          MatSetValues(M_VSTAR_FIXED, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          bc_contribution = VISC * area_f * alpha_12 * (grad_v_star[cell_id1-1].bc_x * n2.x() + grad_v_star[cell_id1-1].bc_y * n2.y());
          bb_vstar[r1] += bc_contribution;
          bb_vstar[r2] -= bc_contribution;
          // (-1.0) MP
          if (grad_p[cell_id1-1].size != size1) {std::cerr<<"size != size1\n"; exit(1);}
          //for (int i = 0; i < size1+1; i++) { val1[i] /= VISC; nval1[i] = -val1[i]; }
          for (int i = 0; i < size1+1; i++)
          {
            col1[i] = grad_p[cell_id1-1].cell_id[i]-1;
            val1[i] = grad_p[cell_id1-1].coef_x[i] * n2.x() + grad_p[cell_id1-1].coef_y[i] * n2.y();
            val1[i] *= (area_f * alpha_12);
            nval1[i] = -val1[i];
          }
          MatSetValues(M_P, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_P, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          bc_contribution = -area_f * alpha_12 * (grad_p[cell_id1-1].bc_x * n2.x() + grad_p[cell_id1-1].bc_y * n2.y());
          bb_p[r1] += bc_contribution;
          bb_p[r2] -= bc_contribution;
          // copy & paste for T (diffusion sign convention matches USTAR/VSTAR, not -MP)
          if (grad_T[cell_id1-1].size != size1) {std::cerr<<"size != size1\n"; exit(1);}
          for (int i = 0; i < size1+1; i++)
          {
            col1[i] = grad_T[cell_id1-1].cell_id[i]-1;
            val1[i] = -grad_T[cell_id1-1].coef_x[i] * n2.x() - grad_T[cell_id1-1].coef_y[i] * n2.y();
            val1[i] *= (DIFF_T * area_f * alpha_12);
            nval1[i] = -val1[i];
          }
          MatSetValues(M_T_FIXED, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_T_FIXED, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          bc_contribution = DIFF_T * area_f * alpha_12 * (grad_T[cell_id1-1].bc_x * n2.x() + grad_T[cell_id1-1].bc_y * n2.y());
          bb_t[r1] += bc_contribution;
          bb_t[r2] -= bc_contribution;

          int size2 = grad_u_star[cell_id2-1].size;
          if (size2 > 3) {std::cerr<<"size2 > 3\n"; exit(1);}
          PetscInt col2[size2+1]; double val2[size2+1], nval2[size2+1];
          for (int i = 0; i < size2+1; i++)
          {
            col2[i] = grad_u_star[cell_id2-1].cell_id[i]-1;
            val2[i] = -grad_u_star[cell_id2-1].coef_x[i] * n2.x() - grad_u_star[cell_id2-1].coef_y[i] * n2.y();
            val2[i] *= (VISC * area_f * alpha_21);
            nval2[i] = -val2[i];
          }
          MatSetValues(M_USTAR_FIXED, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_USTAR_FIXED, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = VISC * area_f * alpha_21 * (grad_u_star[cell_id2-1].bc_x * n2.x() + grad_u_star[cell_id2-1].bc_y * n2.y());
          bb_ustar[r1] += bc_contribution;
          bb_ustar[r2] -= bc_contribution;
          // copy & paste for VSTAR
          if (grad_v_star[cell_id2-1].size != size2) {std::cerr<<"size != size2\n"; exit(1);}
          for (int i = 0; i < size2+1; i++)
          {
            col2[i] = grad_v_star[cell_id2-1].cell_id[i]-1;
            val2[i] = -grad_v_star[cell_id2-1].coef_x[i] * n2.x() - grad_v_star[cell_id2-1].coef_y[i] * n2.y();
            val2[i] *= (VISC * area_f * alpha_21);
            nval2[i] = -val2[i];
          }
          MatSetValues(M_VSTAR_FIXED, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_VSTAR_FIXED, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = VISC * area_f * alpha_21 * (grad_v_star[cell_id2-1].bc_x * n2.x() + grad_v_star[cell_id2-1].bc_y * n2.y());
          bb_vstar[r1] += bc_contribution;
          bb_vstar[r2] -= bc_contribution;
          // (-1.0) MP
          if (grad_p[cell_id2-1].size != size2) {std::cerr<<"size != size2\n"; exit(1);}
          for (int i = 0; i < size2+1; i++)
          {
            col2[i] = grad_p[cell_id2-1].cell_id[i]-1;
            val2[i] = grad_p[cell_id2-1].coef_x[i] * n2.x() + grad_p[cell_id2-1].coef_y[i] * n2.y();
            val2[i] *= (area_f * alpha_21);
            nval2[i] = -val2[i];
          }
          //for (int i = 0; i < size2+1; i++) { val2[i] /= VISC; nval2[i] = -val2[i]; }
          MatSetValues(M_P, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_P, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = -area_f * alpha_21 * (grad_p[cell_id2-1].bc_x * n2.x() + grad_p[cell_id2-1].bc_y * n2.y());
          bb_p[r1] += bc_contribution;
          bb_p[r2] -= bc_contribution;
          // copy & paste for T
          if (grad_T[cell_id2-1].size != size2) {std::cerr<<"size != size2\n"; exit(1);}
          for (int i = 0; i < size2+1; i++)
          {
            col2[i] = grad_T[cell_id2-1].cell_id[i]-1;
            val2[i] = -grad_T[cell_id2-1].coef_x[i] * n2.x() - grad_T[cell_id2-1].coef_y[i] * n2.y();
            val2[i] *= (DIFF_T * area_f * alpha_21);
            nval2[i] = -val2[i];
          }
          MatSetValues(M_T_FIXED, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_T_FIXED, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = DIFF_T * area_f * alpha_21 * (grad_T[cell_id2-1].bc_x * n2.x() + grad_T[cell_id2-1].bc_y * n2.y());
          bb_t[r1] += bc_contribution;
          bb_t[r2] -= bc_contribution;
      }
    }
  }

  // Time derivative term
  for(std::vector<FluentTriCell>::iterator it = cell_set.begin(); it != cell_set.end(); ++it)
  {
    double RHO = 1.0;

    PetscInt row = it->id() - 1;
    PetscInt col[1];
    double val[1];

    col[0] = row; val[0] = RHO * it->volume() / DT;
    MatSetValues(M_USTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
    bb_ustar[row] += uu[row] * RHO * it->volume() / DT;

    MatSetValues(M_VSTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
    bb_vstar[row] += vv[row] * RHO * it->volume() / DT;
    /* This function has been checked; 09/02/2019, 1:38PM */

    MatSetValues(M_T_FIXED, 1, &row, 1, col, val, ADD_VALUES);
    bb_t[row] += tt[row] * RHO * it->volume() / DT;
  }

  MatAssemblyBegin(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_P, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_T_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_T_FIXED, MAT_FINAL_ASSEMBLY);

  if (PRESSURE_ANCHOR_NEEDED)
  {
    // No zone provides a pressure Dirichlet B.C., so the pressure equation (a pure
    // Neumann/diffusion problem) is only defined up to an additive constant. Anchor
    // cell 0's pressure to zero to make the linear system solvable; the matching
    // bb_p[0] = 0.0 in updateMassVeclocities() keeps this row's RHS consistent.
    PetscInt row[1]; row[0] = 0;
    MatZeroRows(M_P, 1, row, 1.0, PETSC_NULL, PETSC_NULL);
  }

  // Seed the reusable M_USTAR/M_VSTAR from the now-finalized fixed matrices (see the
  // comment by their KSPCreate calls above), and set each KSP's operators/options once
  // here rather than on every solve -- M_P's operator/values never change again after
  // this point, and M_USTAR/M_VSTAR keep the same Mat handle for the rest of the run
  // (only their values change, via MatCopy each Picard iteration), so PETSc picks up
  // those value changes automatically on the next KSPSolve without re-calling
  // KSPSetOperators.
  MatDuplicate(M_USTAR_FIXED, MAT_COPY_VALUES, &M_USTAR);
  MatDuplicate(M_VSTAR_FIXED, MAT_COPY_VALUES, &M_VSTAR);
  MatDuplicate(M_T_FIXED, MAT_COPY_VALUES, &M_T);
  KSPSetOperators(ksp_USTAR, M_USTAR, M_USTAR);
  KSPSetFromOptions(ksp_USTAR);
  KSPSetOperators(ksp_VSTAR, M_VSTAR, M_VSTAR);
  KSPSetFromOptions(ksp_VSTAR);
  KSPSetOperators(ksp_P, M_P, M_P);
  KSPSetFromOptions(ksp_P);
  KSPSetOperators(ksp_T, M_T, M_T);
  KSPSetFromOptions(ksp_T);

  VecRestoreArray(b_USTAR, &bb_ustar);
  VecRestoreArray(b_VSTAR, &bb_vstar);
  VecRestoreArray(b_p, &bb_p);
  VecRestoreArray(b_T, &bb_t);
  VecRestoreArray(u, &uu);
  VecRestoreArray(v, &vv);
  VecRestoreArray(T, &tt);
  //VecRestoreArray(u_STAR, &uu_star);
  //VecRestoreArray(v_STAR, &vv_star);
  // std::cout << "KSP START" << std::endl;

  unsigned int N_TIMESTEP = (unsigned int)n_timestep_opt;

  // Monolithic Newton-Krylov path ("-newton true"). Everything above this point is
  // shared setup -- mesh, B.C.s, the GRAD stencils, and the assembled diffusion +
  // time-derivative operators -- all of which the residual in newton_solver.C reuses
  // rather than re-deriving. It runs its own time loop and returns with u/v/p/T at the
  // final time level, so the cleanup below is unchanged.
  if (NEWTON_ON)
  {
    NewtonSolverConfig cfg;
    cfg.p_mesh = app.p_mesh;
    cfg.grad_p = &grad_p;
    cfg.M_USTAR_FIXED = M_USTAR_FIXED; cfg.M_VSTAR_FIXED = M_VSTAR_FIXED; cfg.M_T_FIXED = M_T_FIXED;
    cfg.M_USTAR = M_USTAR; cfg.M_VSTAR = M_VSTAR; cfg.M_T = M_T;
    cfg.b_USTAR = b_USTAR; cfg.b_VSTAR = b_VSTAR; cfg.b_T = b_T;
    cfg.u = u; cfg.v = v; cfg.p = p; cfg.T = T;
    cfg.F_face = F_face_star;
    cfg.n_timestep = (int)N_TIMESTEP;
    cfg.output_interval = (int)output_interval;
    cfg.exoCtx = exoCtx;
    cfg.convergence_file = &convergence_file;

    runNewtonSolver(cfg);

    closeExodusOutput(exoCtx);
    VecDestroy(&u); VecDestroy(&v); VecDestroy(&p); VecDestroy(&u_STAR); VecDestroy(&v_STAR);
    VecDestroy(&b_USTAR); VecDestroy(&b_VSTAR); VecDestroy(&b_p); VecDestroy(&p_old_it);
    VecDestroy(&b_T); VecDestroy(&b_T_FINAL); VecDestroy(&T); VecDestroy(&T_NEW);
    VecDestroy(&b_USTAR_FINAL); VecDestroy(&b_VSTAR_FINAL); VecDestroy(&b_p_FINAL);
    VecDestroy(&grad_p_x); VecDestroy(&grad_p_y); VecDestroy(&p_src_x); VecDestroy(&p_src_y);
    VecDestroy(&F_face_star); VecDestroy(&F_0f_star); VecDestroy(&F_face_star_old_it);
    MatDestroy(&M_USTAR_FIXED); MatDestroy(&M_VSTAR_FIXED); MatDestroy(&M_USTAR); MatDestroy(&M_VSTAR);
    MatDestroy(&M_P); MatDestroy(&M_T_FIXED); MatDestroy(&M_T);
    KSPDestroy(&ksp_USTAR); KSPDestroy(&ksp_VSTAR); KSPDestroy(&ksp_P); KSPDestroy(&ksp_T);
    app.FreeWorkSpace();
    PetscFinalize();
    return 0;
  }

  // Index of the next Exodus time record. ExodusII records must be contiguous (see
  // DMSetOutputSequenceNumber/ex_put_time in ExodusOutput.C), so this counts only the
  // steps actually written and parts company with time_step once output_interval > 1.
  int output_index = 0;
  for (int time_step = 0; time_step < N_TIMESTEP; time_step++)
  {
    // Begin of time step
    // Perform u_star/v_star -> F_0f -> P -> F_f iteration to get a converged F_f (as well as P)
    std::cout << "Time Step " << time_step << " ... ";
    int N_it_max = 100;
    for (int i = 0; i < N_it_max; i++)
    {
      MatCopy(M_USTAR_FIXED, M_USTAR, SAME_NONZERO_PATTERN);
      MatCopy(M_VSTAR_FIXED, M_VSTAR, SAME_NONZERO_PATTERN);

      // update advection operator to solve u_star and v_star
      VecCopy(b_USTAR, b_USTAR_FINAL);
      VecCopy(b_VSTAR, b_VSTAR_FINAL);
      updateAdvectionOperator(app.p_mesh, F_face_star, b_USTAR_FINAL, b_VSTAR_FINAL, M_USTAR, M_VSTAR);
      // Boussinesq buoyancy (S_v = T, lagged one time step): added fresh into
      // b_VSTAR_FINAL every time it's rebuilt, never into the persistent b_VSTAR
      // baseline, so it never accumulates across time steps or inner iterations.
      if (BUOYANCY_ON)
        updateBuoyancySource(app.p_mesh, T, b_VSTAR_FINAL);

      KSPSolve(ksp_USTAR, b_USTAR_FINAL, u_STAR);
      KSPSolve(ksp_VSTAR, b_VSTAR_FINAL, v_STAR);

      // Update "mass velocity" flux F_0f_star, after u_star and v_star are solved
      VecCopy(b_p, b_p_FINAL);
      updateMassVeclocities(app.p_mesh, u_STAR, v_STAR, F_0f_star, b_p_FINAL);

      // Solve pressure equation
      KSPSolve(ksp_P, b_p_FINAL, p);

      // Compute grad_P from solved p
      evaluatePressureGradientValues(app.p_mesh, p, grad_p_x, grad_p_y, grad_p.data());

      // update F_face_star from solved pressure, and check if it has converged
      VecCopy(F_face_star, F_face_star_old_it);
      updateFfaceStar(app.p_mesh, F_face_star, F_0f_star, p, grad_p_x, grad_p_y);
      VecAXPY(F_face_star_old_it, -1.0, F_face_star); // now F_face_star_old_it is F_face_star - F_face_star_old_it
      PetscScalar error_norm = 0;
      VecNorm(F_face_star_old_it, NORM_2, &error_norm);

      if (error_norm < 1.0e-9)
      {
        // std::cout << "it = " << i << ": L2 error = " << error_norm << std::endl;
        break;
      }
      if (i == N_it_max - 1)
      {
        std::cout << "it = " << i << ": L2 error = " << error_norm << ". Inner iterations did not converge." << std::endl;
        exit(1);
      }
    }
    // After F_f and P are converged, solve u_star and v_star again as the final solutions for u and v
    MatCopy(M_USTAR_FIXED, M_USTAR, SAME_NONZERO_PATTERN);
    MatCopy(M_VSTAR_FIXED, M_VSTAR, SAME_NONZERO_PATTERN);
    // update advection operator to solve u_star and v_star
    VecCopy(b_USTAR, b_USTAR_FINAL);
    VecCopy(b_VSTAR, b_VSTAR_FINAL);
    updateAdvectionOperator(app.p_mesh, F_face_star, b_USTAR_FINAL, b_VSTAR_FINAL, M_USTAR, M_VSTAR);
    if (BUOYANCY_ON)
      updateBuoyancySource(app.p_mesh, T, b_VSTAR_FINAL);
    // update rhs due to pressure gradient
    updatePressureGradientAsSource(app.p_mesh, b_USTAR_FINAL, b_VSTAR_FINAL, grad_p_x, grad_p_y);

    // Solve u_star and v_star considering pressure gradients
    KSPSolve(ksp_USTAR, b_USTAR_FINAL, u_STAR);  // now u_star is u^(n+1)
    KSPSolve(ksp_VSTAR, b_VSTAR_FINAL, v_STAR);  // now v_star is v^(n+1)

    // Solve energy equation once per time step, using the now-converged F_face_star
    // this means that the temperature effect is not fed back to y-momentum equation, such that
    // there is a time lag, and that when Boussinesq buoyancy is considered, this becomes a
    // psuedo-transient, instead of true transient
    MatCopy(M_T_FIXED, M_T, SAME_NONZERO_PATTERN);
    VecCopy(b_T, b_T_FINAL);
    updateEnergyAdvectionOperator(app.p_mesh, F_face_star, b_T_FINAL, M_T);
    KSPSolve(ksp_T, b_T_FINAL, T_NEW);  // now T_NEW is T^(n+1)

    // Update rhs of velocity equations, prepare for next time step
    //VecAXPY(b_USTAR, -1.0, p_src_x);
    //VecAXPY(b_VSTAR, -1.0, p_src_y);
    PetscScalar * uu_star, * vv_star;
    VecGetArray(b_USTAR, &bb_ustar);  VecGetArray(b_VSTAR, &bb_vstar);
    VecGetArray(u, &uu);  VecGetArray(v, &vv);  VecGetArray(u_STAR, &uu_star);  VecGetArray(v_STAR, &vv_star);
    for(std::vector<FluentTriCell>::iterator it = cell_set.begin(); it != cell_set.end(); ++it)
    {
      //double DT = 0.01;
      double RHO = 1.0;

      PetscInt row = it->id() - 1;

      bb_ustar[row] += (uu_star[row] - uu[row]) * RHO * it->volume() / DT;
      bb_vstar[row] += (vv_star[row] - vv[row]) * RHO * it->volume() / DT;
    }
    VecRestoreArray(b_USTAR, &bb_ustar);  VecRestoreArray(b_VSTAR, &bb_vstar);
    VecRestoreArray(u, &uu);  VecRestoreArray(v, &vv);  VecRestoreArray(u_STAR, &uu_star);  VecRestoreArray(v_STAR, &vv_star);

    // Update rhs of energy equation, prepare for next time step (mirrors b_USTAR/b_VSTAR above)
    PetscScalar * tt_new;
    VecGetArray(b_T, &bb_t);  VecGetArray(T, &tt);  VecGetArray(T_NEW, &tt_new);
    for(std::vector<FluentTriCell>::iterator it = cell_set.begin(); it != cell_set.end(); ++it)
    {
      double RHO = 1.0;
      PetscInt row = it->id() - 1;
      bb_t[row] += (tt_new[row] - tt[row]) * RHO * it->volume() / DT;
    }
    VecRestoreArray(b_T, &bb_t);  VecRestoreArray(T, &tt);  VecRestoreArray(T_NEW, &tt_new);

    // update u and v from u_star and v_star
    /*check time step convergence first*/
    VecAXPY(u, -1.0, u_STAR); // now u is (u_STAR - u_old)
    VecAXPY(v, -1.0, v_STAR); // now v is (v_STAR - v_old)
    PetscScalar error_norm_u, error_norm_v, error_norm_T;
    VecNorm(u, NORM_2, &error_norm_u); VecNorm(v, NORM_2, &error_norm_v);
    // std::cout << "    u L2 error = " << error_norm_u << std::endl;
    // std::cout << "    v L2 error = " << error_norm_v << std::endl;
    VecAXPY(T, -1.0, T_NEW); // now T is (T_NEW - T_old)
    VecNorm(T, NORM_2, &error_norm_T);
    convergence_file << time_step << "," << error_norm_u << "," << error_norm_v << "," << error_norm_T << std::endl;

    VecCopy(u_STAR, u);
    VecCopy(v_STAR, v);
    VecCopy(T_NEW, T);
    // End of time step

    // The last step is always written, whatever the interval: it carries the stream
    // function and, for the steady-state runs, it is the answer.
    bool computeStreamFunction = (time_step == N_TIMESTEP-1);
    if ((time_step + 1) % output_interval == 0 || computeStreamFunction)
    {
      // drop off VTU support
      // writeOutputFile(output_prefix, output_index, app.p_mesh, u, v, p, F_face_star, computeStreamFunction);
      writeExodusStep(exoCtx, app.p_mesh, output_index, (time_step+1)*DT, u, v, p, T, F_face_star, computeStreamFunction);
      output_index++;
    }
    // std::cout << "Time Step " << time_step << " end <<<<<<<<<" << std::endl << std::endl;
    std::cout << "  [ok]\n";
  }

  closeExodusOutput(exoCtx);

  VecDestroy(&b_USTAR); VecDestroy(&u_STAR); VecDestroy(&u); MatDestroy(&M_USTAR_FIXED); MatDestroy(&M_USTAR);
  VecDestroy(&b_VSTAR); VecDestroy(&v_STAR); VecDestroy(&v); MatDestroy(&M_VSTAR_FIXED); MatDestroy(&M_VSTAR);
  VecDestroy(&p); VecDestroy(&b_p); VecDestroy(&p_old_it); MatDestroy(&M_P);
  VecDestroy(&b_T); VecDestroy(&b_T_FINAL); VecDestroy(&T); VecDestroy(&T_NEW); MatDestroy(&M_T_FIXED); MatDestroy(&M_T);
  KSPDestroy(&ksp_USTAR); KSPDestroy(&ksp_VSTAR); KSPDestroy(&ksp_P); KSPDestroy(&ksp_T);

  VecDestroy(&F_face_star); VecDestroy(&F_0f_star); VecDestroy(&F_face_star_old_it);
  VecDestroy(&grad_p_x); VecDestroy(&grad_p_y); VecDestroy(&p_src_x); VecDestroy(&p_src_y);

  app.FreeWorkSpace();

  PetscFinalize();

  return 0;
}
