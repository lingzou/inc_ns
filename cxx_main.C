#include <iostream>
#include <petscsnes.h>

#include "inc_problems.h"



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

  PETSC_APPCTX app;

  app.p_mesh = new FluentTwoDMesh();
  app.p_mesh->createMeshFromFile("cavity.msh", false, true);
  std::cout << "# of faces: " << app.p_mesh->n_Faces() << std::endl;

  const std::vector<Node> & node_set = app.p_mesh->getNodeSet();
  std::map<int, std::vector<Face> > & face_zone_map = app.p_mesh->getFaceZoneMap();
  std::vector<FluentTriCell> & cell_set = app.p_mesh->getCellSet();

  std::cout << "Number of nodes = " << node_set.size() << std::endl;
  std::cout << "Number of cells = " << cell_set.size() << std::endl;

  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    long int zone_id = it->first;
    std::vector<Face> & faces = it->second;

    std::cout << "Face Zone " << zone_id << std::endl;
    std::cout << "   number of faces = " << faces.size() << std::endl;
  }

  long int n_Cell = app.p_mesh->n_Cells();
  long int n_Face = app.p_mesh->n_Faces();

  /* Setup B.C. */
  std::map<int, double> U_BC;
  std::map<int, double> V_BC;
  U_BC[2] = 0.0; U_BC[3] = 0.0; U_BC[4] = 0.0; U_BC[5] = 1.0;
  V_BC[2] = 0.0; V_BC[3] = 0.0; V_BC[4] = 0.0; V_BC[5] = 0.0;
  /* End of B.C. setting */

  // --->
  std::cout << "GRAD START" << std::endl;
  GRAD grad_u_star[n_Cell];
  GRAD grad_v_star[n_Cell];
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    std::cout << "zone = " << zone << std::endl;
    switch (zone)
    {
      case 2: // RIGHT
      case 3: // BOTTOM
      case 4: // LEFT
      case 5: // TOP
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          Vec3d face_normal = face.faceNormal();
          double v1 = cell_set.at(cell_id1-1).volume();

          grad_u_star[cell_id1-1].addBCContribution(U_BC[zone] * face_normal.x() / v1, U_BC[zone] * face_normal.y() / v1);
          grad_v_star[cell_id1-1].addBCContribution(V_BC[zone] * face_normal.x() / v1, V_BC[zone] * face_normal.y() / v1);
        }
      }
      break;

      case 7:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();
          double v1 = cell_set.at(cell_id1-1).volume();
          double v2 = cell_set.at(cell_id2-1).volume();
          Vec3d face_normal = face.faceNormal();

          double alpha_12 = v2 / (v1 + v2);
          double alpha_21 = 1.0 - alpha_12;
          // Cell 1
          grad_u_star[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);
          grad_v_star[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);

          grad_u_star[cell_id1-1].size++;
          grad_v_star[cell_id1-1].size++;

          grad_u_star[cell_id1-1].addCoef(grad_u_star[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
          grad_v_star[cell_id1-1].addCoef(grad_v_star[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
          // Cell 2
          grad_u_star[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);
          grad_v_star[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);

          grad_u_star[cell_id2-1].size++;
          grad_v_star[cell_id2-1].size++;

          grad_u_star[cell_id2-1].addCoef(grad_u_star[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
          grad_v_star[cell_id2-1].addCoef(grad_v_star[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  std::cout << "GRAD END" << std::endl;
  // <---

  // KSP stuffs
  std::cout << "KSP START" << std::endl;
  Mat       M_USTAR_FIXED, M_VSTAR_FIXED;
  Mat       M_USTAR, M_VSTAR;
  Mat       M_P;
  Vec       u, u_STAR, b_USTAR;
  Vec       v, v_STAR, b_VSTAR;
  Vec       p, b_p, p_old_it;
  Vec       F_face_star, F_0f_star, F_face_star_old_it;
  KSP       ksp_USTAR, ksp_VSTAR, ksp_P;

  Vec       p_src_x, p_src_y;

  std::cout << "KSP START SETUP" << std::endl;

  VecCreate(PETSC_COMM_SELF, &u_STAR);
  VecSetSizes(u_STAR, PETSC_DECIDE, n_Cell);
  VecSetFromOptions(u_STAR);
  VecDuplicate(u_STAR, &b_USTAR);
  VecDuplicate(u_STAR, &v_STAR);
  VecDuplicate(u_STAR, &b_VSTAR);
  VecDuplicate(u_STAR, &u);
  VecDuplicate(u_STAR, &v);
  VecDuplicate(u_STAR, &p);
  VecDuplicate(u_STAR, &b_p);
  VecDuplicate(u_STAR, &p_old_it);

  VecDuplicate(u_STAR, &p_src_x);
  VecDuplicate(u_STAR, &p_src_y);

  VecSet(u_STAR, 0.0);
  VecSet(v_STAR, 0.0);
  VecSet(b_USTAR, 0.0);
  VecSet(b_VSTAR, 0.0);
  VecSet(u, 0.0);
  VecSet(v, 0.0);
  VecSet(p, 0.0);
  VecSet(b_p, 0.0);
  VecSet(p_old_it, 0.0);

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
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_USTAR);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_USTAR, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_VSTAR_FIXED);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_VSTAR_FIXED, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_VSTAR);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_VSTAR, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_P);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_P, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  KSPCreate(PETSC_COMM_WORLD, &ksp_USTAR);
  KSPCreate(PETSC_COMM_WORLD, &ksp_VSTAR);
  KSPCreate(PETSC_COMM_WORLD, &ksp_P);
  //KSPSetType(ksp_P, KSPCG);

  std::cout << "KSP END SETUP" << std::endl;

  PetscScalar * bb_ustar;
  PetscScalar * bb_vstar;
  PetscScalar * uu;
  PetscScalar * vv;
  VecGetArray(b_USTAR, &bb_ustar);
  VecGetArray(b_VSTAR, &bb_vstar);
  VecGetArray(u, &uu);
  VecGetArray(v, &vv);

  // Diffusion terms
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    //double VISC = 1.0;
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        std::cout << "zone = " << zone << std::endl;
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second).at(j);
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          const FluentTriCell & cell = cell_set.at(cell_id-1);

          long int node_id1 = face.node_id1();
          long int node_id2 = face.node_id2();

          Point pt0 = cell.centroid();
          Point pt1 = node_set.at(node_id1-1).point();
          Point pt2 = node_set.at(node_id2-1).point();
          double distance = app.p_mesh->node_to_face_distance(pt0, pt1, pt2);

          PetscInt row = cell_id - 1;
          PetscInt col[1]; col[0] = row;
          double val[1];
          val[0] = VISC * face.area() / distance;
          MatSetValues(M_USTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
          MatSetValues(M_VSTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);

          bb_ustar[row] += VISC * face.area() / distance * U_BC[zone];
          bb_vstar[row] += VISC * face.area() / distance * V_BC[zone];
        }
      }
      break;

      case 7:
      {
        std::cout << "zone = " << zone << std::endl;
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second).at(j);
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          const FluentTriCell & cell_1 = cell_set.at(cell_id1-1);
          const FluentTriCell & cell_2 = cell_set.at(cell_id2-1);
          double v1 = cell_1.volume();
          double v2 = cell_2.volume();
          const Point & ct1 = cell_1.centroid();
          const Point & ct2 = cell_2.centroid();
          Vec3d ct_to_ct = ct2 - ct1;
          double distance = ct_to_ct.norm();
          Vec3d n1 = ct_to_ct.unitVector();

          Vec3d face_normal = face.faceNormal();
          Vec3d nf = face_normal.unitVector();
          Vec3d n2 = nf - n1;

          double area_f = face.area();
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

          col[0] = r2; val[0] = VISC * area_f / distance; nval[0] = -val[0];
          MatSetValues(M_USTAR_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag
          // copy & paste for VSTAR and (-1.0)MP
          MatSetValues(M_VSTAR_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag
          val[0] = area_f / distance; nval[0] = -val[0];
          MatSetValues(M_P, 1, &r1, 1, col, val, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_P, 1, &r2, 1, col, nval, ADD_VALUES);   // cell 2 diag

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
          MatSetValues(M_VSTAR_FIXED, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          bc_contribution = VISC * area_f * alpha_12 * (grad_v_star[cell_id1-1].bc_x * n2.x() + grad_v_star[cell_id1-1].bc_y * n2.y());
          bb_vstar[r1] += bc_contribution;
          bb_vstar[r2] -= bc_contribution;
          // (-1.0) MP
          for (int i = 0; i < size1+1; i++) { val1[i] /= VISC; nval1[i] = -val1[i]; }
          MatSetValues(M_P, 1, &r1, size1+1, col1, nval1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_P, 1, &r2, size1+1, col1, val1, ADD_VALUES); // cell 2 off-diag

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
          MatSetValues(M_VSTAR_FIXED, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_VSTAR_FIXED, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = VISC * area_f * alpha_21 * (grad_v_star[cell_id2-1].bc_x * n2.x() + grad_v_star[cell_id2-1].bc_y * n2.y());
          bb_vstar[r1] += bc_contribution;
          bb_vstar[r2] -= bc_contribution;
          // (-1.0) MP
          for (int i = 0; i < size2+1; i++) { val2[i] /= VISC; nval2[i] = -val2[i]; }
          MatSetValues(M_P, 1, &r1, size2+1, col2, nval2, ADD_VALUES);
          MatSetValues(M_P, 1, &r2, size2+1, col2, val2, ADD_VALUES);
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
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
  }

  MatAssemblyBegin(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_P, MAT_FINAL_ASSEMBLY);
  PetscInt row[1]; row[0] = 0; // anchor p[0] = 0.0
  MatZeroRows(M_P, 1, row, 1.0, PETSC_NULL, PETSC_NULL);

  VecRestoreArray(b_USTAR, &bb_ustar);
  VecRestoreArray(b_VSTAR, &bb_vstar);
  VecRestoreArray(u, &uu);
  VecRestoreArray(v, &vv);
  // std::cout << "KSP START" << std::endl;

  for (int time_step = 0; time_step < 20; time_step++)
  {
    // Begin of time step
    // Perform u_star/v_star -> F_0f -> P -> F_f iteration to get a converged F_f (as well as P)
    for (int i = 0; i < 20; i++)
    {
      MatDuplicate(M_USTAR_FIXED, MAT_COPY_VALUES, &M_USTAR);
      MatDuplicate(M_VSTAR_FIXED, MAT_COPY_VALUES, &M_VSTAR);

      // update advection operator to solve u_star and v_star
      updateAdvectionOperator(app.p_mesh, F_face_star, M_USTAR, M_VSTAR);

      KSPSetOperators(ksp_USTAR, M_USTAR, M_USTAR);
      KSPSetFromOptions(ksp_USTAR);
      KSPSetOperators(ksp_VSTAR, M_VSTAR, M_VSTAR);
      KSPSetFromOptions(ksp_VSTAR);
      KSPSolve(ksp_USTAR, b_USTAR, u_STAR);
      KSPSolve(ksp_VSTAR, b_VSTAR, v_STAR);

      // Update "mass velocity" flux F_0f_star, after u_star and v_star are solved
      updateMassVeclocities(app.p_mesh, u_STAR, v_STAR, F_0f_star, b_p);

      // Solve pressure equation
      KSPSetOperators(ksp_P, M_P, M_P);
      KSPSetFromOptions(ksp_P);
      //VecCopy(p, p_old_it);
      KSPSolve(ksp_P, b_p, p);
      /*
      VecAXPY(p_old_it, -1.0, p); // now p_old_it is p - p_old_it
      PetscScalar error_norm, error_inf;
      VecNorm(p_old_it, NORM_2, &error_norm);
      VecNorm(p_old_it, NORM_INFINITY, &error_inf);
      std::cout << "p error_norm = " << error_norm << std::endl;
      std::cout << "p error_inf = " << error_inf << std::endl;*/

      // update F_face_star from solved pressure
      VecCopy(F_face_star, F_face_star_old_it);
      updateFfaceStar(app.p_mesh, F_face_star, F_0f_star, p, grad_u_star);
      VecAXPY(F_face_star_old_it, -1.0, F_face_star); // now F_face_star_old_it is F_face_star - F_face_star_old_it
      PetscScalar error_norm;
      VecNorm(F_face_star_old_it, NORM_2, &error_norm);
      //VecNorm(F_face_star_old_it, NORM_INFINITY, &error_inf);
      //std::cout << "F error_norm = " << error_norm << std::endl;
      //std::cout << "F error_inf = " << error_inf << std::endl;

      if (error_norm < 1.0e-9)
      {
        std::cout << "it = " << i << ": L2 error = " << error_norm << std::endl;
        break;
      }
    }
    // After F_f and P are converged, solve u_star and v_star again as the final solutions for u and v
    MatDuplicate(M_USTAR_FIXED, MAT_COPY_VALUES, &M_USTAR);
    MatDuplicate(M_VSTAR_FIXED, MAT_COPY_VALUES, &M_VSTAR);
    // update advection operator to solve u_star and v_star
    updateAdvectionOperator(app.p_mesh, F_face_star, M_USTAR, M_VSTAR);
    // update rhs due to pressure gradient
    updatePressureGradientAsSource(app.p_mesh, p, p_src_x, p_src_y);
    VecAXPY(b_USTAR, 1.0, p_src_x);
    VecAXPY(b_VSTAR, 1.0, p_src_y);

    KSPSetOperators(ksp_USTAR, M_USTAR, M_USTAR);
    KSPSetFromOptions(ksp_USTAR);
    KSPSetOperators(ksp_VSTAR, M_VSTAR, M_VSTAR);
    KSPSetFromOptions(ksp_VSTAR);
    KSPSolve(ksp_USTAR, b_USTAR, u_STAR);  // now u_star is u^(n+1)
    KSPSolve(ksp_VSTAR, b_VSTAR, v_STAR);  // now v_star is v^(n+1)

    // Update rhs of velocity equations, prepare for next time step
    VecAXPY(b_USTAR, -1.0, p_src_x);
    VecAXPY(b_VSTAR, -1.0, p_src_y);
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

    // update u and v from u_star and v_star
    /*check time step convergence first*/
    VecAXPY(u, -1.0, u_STAR); // now u is (u_STAR - u_old)
    VecAXPY(v, -1.0, v_STAR); // now v is (v_STAR - v_old)
    PetscScalar error_norm_u, error_norm_v;
    VecNorm(u, NORM_2, &error_norm_u); VecNorm(v, NORM_2, &error_norm_v);
    std::cout << "    u L2 error = " << error_norm_u << std::endl;
    std::cout << "    v L2 error = " << error_norm_v << std::endl;

    VecCopy(u_STAR, u);
    VecCopy(v_STAR, v);
    // End of time step

    std::string file_name = "output/output_" + std::to_string(time_step) + ".vtu";
    FILE * ptr_File;
    ptr_File = fopen(file_name.c_str(), "w");
    app.p_mesh->writeMesh(ptr_File);
    std::ostringstream out_string_stream;
    out_string_stream << "      <CellData>" << "\n";

    // CELL DATA (cell ID)
    out_string_stream << "        <DataArray type=\"Float32\" Name=\"Cell_ID\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < cell_set.size(); i++)
      out_string_stream << "          " << cell_set[i].id() << "\n";
    out_string_stream << "        </DataArray>" << "\n";
    // CELL DATA (volume)
    out_string_stream << "        <DataArray type=\"Float32\" Name=\"volume\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < cell_set.size(); i++)
      out_string_stream << "          " << cell_set[i].volume() << "\n";
    out_string_stream << "        </DataArray>" << "\n";

    //PetscScalar * uu;
    VecGetArray(u, &uu);
    out_string_stream << "        <DataArray type=\"Float32\" Name=\"u_star\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < cell_set.size(); i++)
      out_string_stream << "          " << uu[i] << "\n";
    out_string_stream << "        </DataArray>" << "\n";
    VecRestoreArray(u, &uu);

    //PetscScalar * vv;
    VecGetArray(v, &vv);
    out_string_stream << "        <DataArray type=\"Float32\" Name=\"v_star\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < cell_set.size(); i++)
      out_string_stream << "          " << vv[i] << "\n";
    out_string_stream << "        </DataArray>" << "\n";
    VecRestoreArray(v, &vv);

    PetscScalar * pp;
    VecGetArray(p, &pp);
    out_string_stream << "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < cell_set.size(); i++)
      out_string_stream << "          " << pp[i] << "\n";
    out_string_stream << "        </DataArray>" << "\n";
    VecRestoreArray(p, &pp);

    out_string_stream << "      </CellData>" << "\n";

    // POINT DATA
    out_string_stream << "      <PointData>" << "\n";
    // NODE ID
    out_string_stream << "        <DataArray type=\"Float32\" Name=\"Node_ID\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < node_set.size(); i++)
      out_string_stream << "          " << node_set[i].id() << "\n";
    out_string_stream << "        </DataArray>" << "\n";
    out_string_stream << "      </PointData>" << "\n";

    fprintf(ptr_File, "%s", out_string_stream.str().c_str());

    app.p_mesh->finishFile(ptr_File);
    fclose(ptr_File);
    //delete p_mesh;
  }


  VecDestroy(&b_USTAR); VecDestroy(&u_STAR); VecDestroy(&u); MatDestroy(&M_USTAR_FIXED); MatDestroy(&M_USTAR);
  VecDestroy(&b_VSTAR); VecDestroy(&v_STAR); VecDestroy(&v); MatDestroy(&M_VSTAR_FIXED); MatDestroy(&M_VSTAR);
  VecDestroy(&p); VecDestroy(&b_p); VecDestroy(&p_old_it); MatDestroy(&M_P);
  KSPDestroy(&ksp_USTAR); KSPDestroy(&ksp_VSTAR); KSPDestroy(&ksp_P);

  VecDestroy(&F_face_star); VecDestroy(&F_0f_star); VecDestroy(&F_face_star_old_it);
  VecDestroy(&p_src_x); VecDestroy(&p_src_y);

  app.FreeWorkSpace();

  PetscFinalize();

  return 0;
}
