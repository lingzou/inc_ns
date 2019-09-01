#include <iostream>
#include <petscsnes.h>

#include "FluentTwoDMesh.h"

PetscErrorCode SNESFormFunction(SNES, Vec, Vec, void *);

struct GRAD
{
  int size;
  long int cell_id[4];
  double coef_x[4], coef_y[4];
  double bc_x, bc_y;
  GRAD()
  {
    size = 0;
    cell_id[0] = -1; cell_id[1] = -1; cell_id[2] = -1; cell_id[3] = -1;
    coef_x[0] = 0.0; coef_x[1] = 0.0; coef_x[2] = 0.0; coef_x[3] = 0.0;
    coef_y[0] = 0.0; coef_y[1] = 0.0; coef_y[2] = 0.0; coef_y[3] = 0.0;
    bc_x = 0.0; bc_y = 0.0;
  }
  void addCoef(int loc, long int id, double x_val, double y_val)
  {
    cell_id[loc] = id;
    coef_x[loc] += x_val;
    coef_y[loc] += y_val;
  }
  void addBCContribution(double x_val, double y_val) { bc_x += x_val; bc_y += y_val; }
};

struct PETSC_APPCTX
{
  Vec T_cell, res;
  Vec T_face;
  Vec grad_T_x, grad_T_y;
  FluentTwoDMesh * p_mesh;
};

PetscErrorCode SNESFormFunction(SNES snes, Vec u, Vec f, void * appCtx)
{
  PETSC_APPCTX * p_appCtx = (PETSC_APPCTX *) appCtx;
  FluentTwoDMesh * p_mesh = p_appCtx->p_mesh;
  std::map<int, std::vector<Face> > & face_zone_map = p_mesh->getFaceZoneMap();
  const std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  const std::vector<Node> & node_set = p_mesh->getNodeSet();

  PetscScalar * T_cc, * T_ff, * gradTx, * gradTy, * res;
  VecGetArray(u, &T_cc);
  VecGetArray(p_appCtx->grad_T_x, &gradTx);
  VecGetArray(p_appCtx->grad_T_y, &gradTy);
  VecGetArray(p_appCtx->T_face, &T_ff);
  VecGetArray(f, &res);

  for (int i = 0; i < cell_set.size(); i++)
  { gradTx[i] = 0.0; gradTy[i] = 0.0;

    const FluentTriCell & cell = cell_set.at(i);
    double xx = cell.centroid().x();
    double yy = cell.centroid().y();
    double pi = 3.1415926;
    double qq = 2.0 * pi * pi * sin(pi * xx) * sin(pi * yy);
    res[i] = -qq * cell.volume();
  }
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        for (int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second).at(j);
          long int face_id = face.id();
          T_ff[face_id] = 0.0;
        }
      }
      break;

      case 7:
      {
        for (int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second).at(j);
          long int face_id = face.id();
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          const FluentTriCell & cell_1 = cell_set.at(cell_id1-1);
          const FluentTriCell & cell_2 = cell_set.at(cell_id2-1);
          double v1 = cell_1.volume();
          double v2 = cell_2.volume();
          Vec3d face_normal = face.faceNormal();

          double T_c1 = T_cc[cell_id1-1];
          double T_c2 = T_cc[cell_id2-1];
          double T_f = (T_c1 * v2 + T_c2 * v1) / (v1 + v2);
          T_ff[face_id] = T_f;

          // cell 1
          gradTx[cell_id1-1] += T_f * face_normal.x() / v1;
          gradTy[cell_id1-1] += T_f * face_normal.y() / v1;
          // cell 2
          gradTx[cell_id2-1] -= T_f * face_normal.x() / v2;
          gradTy[cell_id2-1] -= T_f * face_normal.y() / v2;
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          const FluentTriCell & cell = cell_set.at(cell_id-1);
          double distance = 1.0 - cell.centroid().y();
          double dTdy = (0.0 - T_cc[cell_id-1]) / distance;
          res[cell_id-1] -= dTdy * face.area();
        }
      }
      break;
      case 2:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          const FluentTriCell & cell = cell_set.at(cell_id-1);
          double distance = 1.0 - cell.centroid().x();
          double dTdx = (0.0 - T_cc[cell_id-1]) / distance;
          res[cell_id-1] -= dTdx * face.area();
        }
      }
      break;
      case 3:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          const FluentTriCell & cell = cell_set.at(cell_id-1);
          double distance = cell.centroid().y() - 0.0;
          double dTdy = (T_cc[cell_id-1] - 0.0) / distance;
          res[cell_id-1] += dTdy * face.area();
        }
      }
      break;
      case 4:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          const FluentTriCell & cell = cell_set.at(cell_id-1);
          double distance = cell.centroid().x() - 0.0;
          double dTdx = (T_cc[cell_id-1] - 0.0) / distance;
          res[cell_id-1] += dTdx * face.area();
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

          double T_c1 = T_cc[cell_id1-1];
          double T_c2 = T_cc[cell_id2-1];
          Vec3d gradT_c1(gradTx[cell_id1-1], gradTy[cell_id1-1], 0.0);
          Vec3d gradT_c2(gradTx[cell_id2-1], gradTy[cell_id2-1], 0.0);

          // diffusion flux
          double dT_dL = (T_c2 - T_c1) / distance;
          double cross_dif = (gradT_c1 * n2) * v2 + (gradT_c2 * n2) * v1;
          cross_dif = cross_dif / (v1 + v2);
          double diff_flux = -(dT_dL + cross_dif) * face.area();

          res[cell_id1-1] += diff_flux;
          res[cell_id2-1] -= diff_flux;
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  VecRestoreArray(u, &T_cc);
  VecRestoreArray(p_appCtx->grad_T_x, &gradTx);
  VecRestoreArray(p_appCtx->grad_T_y, &gradTy);
  VecRestoreArray(p_appCtx->T_face, &T_ff);
  VecRestoreArray(f, &res);

  return 0;
};

int main(int argc, char **argv)
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  FluentTwoDMesh * p_mesh = new FluentTwoDMesh();
  p_mesh->createMeshFromFile("cavity.msh", false, true);
  std::cout << "# of faces: " << p_mesh->n_Faces() << std::endl;

  const std::vector<Node> & node_set = p_mesh->getNodeSet();
  std::map<int, std::vector<Face> > & face_zone_map = p_mesh->getFaceZoneMap();
  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();

  std::cout << "Number of nodes = " << node_set.size() << std::endl;
  std::cout << "Number of cells = " << cell_set.size() << std::endl;

  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    long int zone_id = it->first;
    std::vector<Face> & faces = it->second;

    std::cout << "Face Zone " << zone_id << std::endl;
    std::cout << "   number of faces = " << faces.size() << std::endl;
  }

  long int n_Cell = p_mesh->n_Cells();
  long int n_Face = p_mesh->n_Faces();


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
      case 5: // TOP
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          //long int cell_id2 = face.cell_id2();
          Vec3d face_normal = face.faceNormal();
          double v1 = cell_set.at(cell_id1-1).volume();

          //std::cout << "cell_ids = " << cell_id1 << " " << cell_id2 << std::endl;
          //std::cout << "node_ids = " << face.node_id1() << " " << face.node_id2() << std::endl;
          //face_normal.print();

          //long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          double U_BC = 1.0;
          grad_u_star[cell_id1-1].addBCContribution(U_BC * face_normal.x() / v1, U_BC * face_normal.y() / v1);
          // V_BC = 0.0, no need to do anything
        }
      }
      break;

      case 2: // RIGHT
      case 3: // BOTTOM
      case 4: // LEFT
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
  Vec       p, b_p;
  Vec       F_face_star, F_0f_star;
  KSP       ksp_USTAR, ksp_VSTAR, ksp_P;

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

  VecSet(u_STAR, 0.0);
  VecSet(v_STAR, 0.0);
  VecSet(b_USTAR, 0.0);
  VecSet(b_VSTAR, 0.0);
  VecSet(u, 0.0);
  VecSet(v, 0.0);
  VecSet(p, 0.0);
  VecSet(b_p, 0.0);

  VecCreate(PETSC_COMM_SELF, &F_face_star);
  VecSetSizes(F_face_star, PETSC_DECIDE, n_Face);
  VecSetFromOptions(F_face_star);
  VecDuplicate(F_face_star, &F_0f_star);
  VecSet(F_face_star, 0.0);
  VecSet(F_0f_star, 0.0);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_USTAR_FIXED);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_USTAR_FIXED, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  KSPCreate(PETSC_COMM_WORLD, &ksp_USTAR);
  //KSPSetOperators(ksp_USTAR, M_USTAR, M_USTAR);
  //KSPSetFromOptions(ksp_USTAR);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_VSTAR_FIXED);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_VSTAR_FIXED, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  KSPCreate(PETSC_COMM_WORLD, &ksp_VSTAR);
  //KSPSetOperators(ksp_VSTAR, M_VSTAR, M_VSTAR);
  //KSPSetFromOptions(ksp_VSTAR);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_Cell, n_Cell, 10, NULL, &M_P);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_P, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  KSPCreate(PETSC_COMM_WORLD, &ksp_P);

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
          double distance = 0.0;

          if (zone == 5)        distance = 1.0 - cell.centroid().y();   // top
          else if (zone == 2)   distance = 1.0 - cell.centroid().x();   // right
          else if (zone == 3)   distance = cell.centroid().y();         // bottom
          else                  distance = cell.centroid().x();         // left

          PetscInt row = cell_id - 1;
          PetscInt col[1]; col[0] = row;
          double val[1];
          val[0] = face.area() / distance;
          MatSetValues(M_USTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
          MatSetValues(M_VSTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);

          double U_BC = 1.0;
          if (zone == 5) bb_ustar[row] = face.area() / distance * U_BC;
          // no rhs for v_star
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

          col[0] = r1; val[0] = area_f / distance; nval[0] = -val[0];
          MatSetValues(M_USTAR_FIXED, 1, &r1, 1, col, val, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, 1, col, nval, ADD_VALUES);  // cell 2 off-diag
          // copy & paste for VSTAR
          MatSetValues(M_VSTAR_FIXED, 1, &r1, 1, col, val, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, 1, col, nval, ADD_VALUES);  // cell 2 off-diag
          col[0] = r2;
          MatSetValues(M_USTAR_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag
          // copy & paste for VSTAR
          MatSetValues(M_VSTAR_FIXED, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag

          int size1 = grad_u_star[cell_id1-1].size;
          if (size1 > 3) {std::cerr<<"size1 > 3\n"; exit(1);}
          PetscInt col1[size1+1]; double val1[size1+1], nval1[size1+1];
          for (int i = 0; i < size1+1; i++)
          {
            col1[i] = grad_u_star[cell_id1-1].cell_id[i]-1;
            val1[i] = -grad_u_star[cell_id1-1].coef_x[i] * n2.x() - grad_u_star[cell_id1-1].coef_y[i] * n2.y();
            val1[i] *= (area_f * alpha_12);
            nval1[i] = -val1[i];
          }
          MatSetValues(M_USTAR_FIXED, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_USTAR_FIXED, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          double bc_contribution = area_f * alpha_12 * (grad_u_star[cell_id1-1].bc_x * n2.x() + grad_u_star[cell_id1-1].bc_y * n2.y());
          bb_ustar[r1] += bc_contribution;
          bb_ustar[r2] -= bc_contribution;
          // copy & paste for VSTAR
          MatSetValues(M_VSTAR_FIXED, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_VSTAR_FIXED, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          bc_contribution = area_f * alpha_12 * (grad_v_star[cell_id1-1].bc_x * n2.x() + grad_v_star[cell_id1-1].bc_y * n2.y());
          bb_vstar[r1] += bc_contribution;
          bb_vstar[r2] -= bc_contribution;

          int size2 = grad_u_star[cell_id2-1].size;
          if (size2 > 3) {std::cerr<<"size2 > 3\n"; exit(1);}
          PetscInt col2[size2+1]; double val2[size2+1], nval2[size2+1];
          for (int i = 0; i < size2+1; i++)
          {
            col2[i] = grad_u_star[cell_id2-1].cell_id[i]-1;
            val2[i] = -grad_u_star[cell_id2-1].coef_x[i] * n2.x() - grad_u_star[cell_id2-1].coef_y[i] * n2.y();
            val2[i] *= (area_f * alpha_21);
            nval2[i] = -val2[i];
          }
          MatSetValues(M_USTAR_FIXED, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_USTAR_FIXED, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = area_f * alpha_21 * (grad_u_star[cell_id2-1].bc_x * n2.x() + grad_u_star[cell_id2-1].bc_y * n2.y());
          bb_ustar[r1] += bc_contribution;
          bb_ustar[r2] -= bc_contribution;
          // copy & paste for VSTAR
          MatSetValues(M_VSTAR_FIXED, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_VSTAR_FIXED, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
          bc_contribution = area_f * alpha_21 * (grad_v_star[cell_id2-1].bc_x * n2.x() + grad_v_star[cell_id2-1].bc_y * n2.y());
          bb_vstar[r1] += bc_contribution;
          bb_vstar[r2] -= bc_contribution;
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }

  MatAssemblyBegin(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatDuplicate(M_VSTAR_FIXED, MAT_COPY_VALUES, &M_P);
  MatScale(M_P, -1.0);

  // Time derivative term
  for(std::vector<FluentTriCell>::iterator it = cell_set.begin(); it != cell_set.end(); ++it)
  {
    double DT = 0.01;
    double RHO = 1.0;

    PetscInt row = it->id() - 1;
    PetscInt col[1];
    double val[1];

    col[0] = row; val[0] = RHO * it->volume() / DT;
    MatSetValues(M_USTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
    bb_ustar[row] += uu[row] * RHO * it->volume() / DT;

    MatSetValues(M_VSTAR_FIXED, 1, &row, 1, col, val, ADD_VALUES);
    bb_vstar[row] += vv[row] * RHO * it->volume() / DT;
  }

  MatAssemblyBegin(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_USTAR_FIXED, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_VSTAR_FIXED, MAT_FINAL_ASSEMBLY);

  VecRestoreArray(b_USTAR, &bb_ustar);
  VecRestoreArray(b_VSTAR, &bb_vstar);
  VecRestoreArray(u, &uu);
  VecRestoreArray(v, &vv);
  std::cout << "KSP START" << std::endl;

for (int i = 0; i < 5; i++)
{
  MatDuplicate(M_USTAR_FIXED, MAT_COPY_VALUES, &M_USTAR);
  MatDuplicate(M_VSTAR_FIXED, MAT_COPY_VALUES, &M_VSTAR);

  //KSPCreate(PETSC_COMM_WORLD, &ksp_USTAR);
  KSPSetOperators(ksp_USTAR, M_USTAR, M_USTAR);
  KSPSetFromOptions(ksp_USTAR);
  //KSPCreate(PETSC_COMM_WORLD, &ksp_VSTAR);
  KSPSetOperators(ksp_VSTAR, M_VSTAR, M_VSTAR);
  KSPSetFromOptions(ksp_VSTAR);
  //KSPCreate(PETSC_COMM_WORLD, &ksp_P);
  KSPSetOperators(ksp_P, M_P, M_P);
  KSPSetFromOptions(ksp_P);

  // Advection terms
  PetscScalar * ff;
  VecGetArray(F_face_star, &ff);
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        // Nothing to do, F_face on all boundaries are zero
      }
      break;

      case 7:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          double GAMMA = 0.5;

          Face & face = (it->second).at(j);
          long int face_id = face.id();
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          const FluentTriCell & cell_1 = cell_set.at(cell_id1-1);
          const FluentTriCell & cell_2 = cell_set.at(cell_id2-1);
          double v1 = cell_1.volume();
          double v2 = cell_2.volume();
          //const Point & ct1 = cell_1.centroid();
          //const Point & ct2 = cell_2.centroid();
          //Vec3d ct_to_ct = ct2 - ct1;
          //double distance = ct_to_ct.norm();
          //Vec3d n1 = ct_to_ct.unitVector();

          //Vec3d face_normal = face.faceNormal();
          //Vec3d nf = face_normal.unitVector();
          //Vec3d n2 = nf - n1;

          double area_f = face.area();
          double alpha_12 = v2 / (v1 + v2);
          double alpha_21 = 1.0 - alpha_12;
          // cell 1
          PetscInt r1 = cell_id1 - 1; PetscInt r2 = cell_id2 - 1;
          PetscInt col[2];
          double val[2], nval[2];

          col[0] = r1; col[1] = r2;
          val[0] = (1.0 - GAMMA) * ff[face_id] * alpha_12;
          val[1] = (1.0 - GAMMA) * ff[face_id] * alpha_21;
          if (ff[face_id] > 0.0)
            val[0] += GAMMA * ff[face_id];
          else
            val[1] += GAMMA * ff[face_id];

          nval[0] = -val[0]; nval[1] = -val[1];
          MatSetValues(M_USTAR, 1, &r1, 2, col, val, ADD_VALUES);
          // copy & paste for VSTAR
          MatSetValues(M_VSTAR, 1, &r1, 2, col, val, ADD_VALUES);

          MatSetValues(M_USTAR, 1, &r2, 2, col, nval, ADD_VALUES);
          // copy & paste for VSTAR
          MatSetValues(M_VSTAR, 1, &r2, 2, col, nval, ADD_VALUES);
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  VecRestoreArray(F_face_star, &ff);
  MatAssemblyBegin(M_USTAR, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_USTAR, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(M_VSTAR, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_VSTAR, MAT_FINAL_ASSEMBLY);

  KSPSolve(ksp_USTAR, b_USTAR, u_STAR);
  KSPSolve(ksp_VSTAR, b_VSTAR, v_STAR);

  // Update "mass velocity" flux F_0f_star
  PetscScalar * f0f;
  PetscScalar * uu_star;
  PetscScalar * vv_star;
  PetscScalar * bb_p;
  VecGetArray(F_0f_star, &f0f);
  VecGetArray(u_STAR, &uu_star);
  VecGetArray(v_STAR, &vv_star);
  VecGetArray(b_p, &bb_p);
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        // Nothing to do, F_face on all boundaries are zero
      }
      break;

      case 7:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second).at(j);
          double r = face.distance_ratio();

          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          double u_face = uu_star[cell_id1-1] * (1.0 - r) + uu_star[cell_id2-1] * r;
          double v_face = vv_star[cell_id1-1] * (1.0 - r) + vv_star[cell_id2-1] * r;

          Vec3d face_normal = face.faceNormal();

          double val = u_face * face_normal.x() + v_face * face_normal.y();
          f0f[face.id()] = val;
          bb_p[cell_id1-1] += val;
          bb_p[cell_id2-1] -= val;

          //std::cout << "j = " << j << " " << face.distance_ratio() << std::endl;
          //std::cout << "j = " << j << " " << face.area() << std::endl;
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }

    // anchor point
    bb_p[0] = 0.0;
  }
  VecRestoreArray(F_0f_star, &f0f);
  VecRestoreArray(u_STAR, &uu_star);
  VecRestoreArray(v_STAR, &vv_star);
  VecRestoreArray(b_p, &bb_p);

  PetscInt row[1];
  row[0] = 0;
  MatZeroRows(M_P, 1, row, 1.0, PETSC_NULL, PETSC_NULL);
  KSPSolve(ksp_P, b_p, p);

  std::cout << "update F_face_star" << std::endl;
  // update F_face_star
  PetscScalar * pp;
  VecGetArray(F_face_star, &ff);
  VecGetArray(F_0f_star, &f0f);
  VecGetArray(p, &pp);
  std::cout << "Start loop" << std::endl;
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        // Nothing to do, F_face on all boundaries are zero
      }
      break;

      case 7:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          //std::cout << "Face j = " << j << std::endl;
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

          //std::cout << "->Cell 1" << std::endl;
          double grad_p_dot_n2_cell1 = 0.0;
          int size1 = grad_u_star[cell_id1-1].size;
          if (size1 > 3) {std::cerr<<"size1 > 3\n"; exit(1);}
          for (int i = 0; i < size1+1; i++)
          {
            long int cell_id_i = grad_u_star[cell_id1-1].cell_id[i]-1;
            double pi = pp[i];
            grad_p_dot_n2_cell1 += (grad_u_star[cell_id1-1].coef_x[i] * n2.x() + grad_u_star[cell_id1-1].coef_y[i] * n2.y()) * pi;
          }
          grad_p_dot_n2_cell1 *= alpha_12;

          //std::cout << "->Cell 2" << std::endl;
          double grad_p_dot_n2_cell2 = 0.0;
          int size2 = grad_u_star[cell_id2-1].size;
          if (size2 > 3) {std::cerr<<"size2 > 3\n"; exit(1);}
          for (int i = 0; i < size2+1; i++)
          {
            long int cell_id_i = grad_u_star[cell_id1-2].cell_id[i]-1;
            double pi = pp[i];
            grad_p_dot_n2_cell2 += (grad_u_star[cell_id2-1].coef_x[i] * n2.x() + grad_u_star[cell_id2-1].coef_y[i] * n2.y()) * pi;
          }
          grad_p_dot_n2_cell2 *= alpha_21;

          double pressure_correction = (pp[cell_id2-1] - pp[cell_id1-1]) / distance + grad_p_dot_n2_cell1 + grad_p_dot_n2_cell2;
          pressure_correction *= area_f;

          double DT = 0.01;
          ff[face.id()] = f0f[face.id()] - DT * pressure_correction;
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  std::cout << "End loop: before restore vec" << std::endl;
  VecRestoreArray(F_face_star, &ff);
  VecRestoreArray(F_0f_star, &f0f);
  VecRestoreArray(p, &pp);
  std::cout << "End loop" << std::endl;
}
  FILE * ptr_File;
  ptr_File = fopen("output/U_STAR.vtu", "w");
  p_mesh->writeMesh(ptr_File);
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

  PetscScalar * uu_star;
  VecGetArray(u_STAR, &uu_star);
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"u_star\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << uu_star[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  VecRestoreArray(u_STAR, &uu_star);

  PetscScalar * vv_star;
  VecGetArray(v_STAR, &vv_star);
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"v_star\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << vv_star[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  VecRestoreArray(v_STAR, &vv_star);

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

  p_mesh->finishFile(ptr_File);
  fclose(ptr_File);
  delete p_mesh;


  VecDestroy(&b_USTAR); VecDestroy(&u_STAR); VecDestroy(&u); MatDestroy(&M_USTAR_FIXED); MatDestroy(&M_USTAR);
  VecDestroy(&b_VSTAR); VecDestroy(&v_STAR); VecDestroy(&v); MatDestroy(&M_VSTAR_FIXED); MatDestroy(&M_VSTAR);
  VecDestroy(&p); VecDestroy(&b_p); MatDestroy(&M_P);
  KSPDestroy(&ksp_USTAR); KSPDestroy(&ksp_VSTAR); KSPDestroy(&ksp_P);

  VecDestroy(&F_face_star); VecDestroy(&F_0f_star);

  PetscFinalize();

  return 0;
}

/*
int main(int argc, char **argv)
{
  // PETSc application starts with PetscInitialize
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  FluentTwoDMesh * p_mesh = new FluentTwoDMesh();
  p_mesh->createMeshFromFile("cavity.msh", false, false);
  std::cout << "# of faces: " << p_mesh->n_Faces() << std::endl;

  const std::vector<Node> & node_set = p_mesh->getNodeSet();
  std::map<int, std::vector<Face> > & face_zone_map = p_mesh->getFaceZoneMap();
  const std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();

  std::cout << "Number of nodes = " << node_set.size() << std::endl;
  std::cout << "Number of cells = " << cell_set.size() << std::endl;

  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    long int zone_id = it->first;
    std::vector<Face> & faces = it->second;

    std::cout << "Face Zone " << zone_id << std::endl;
    std::cout << "   number of faces = " << faces.size() << std::endl;
  }

  long int n_Cell = p_mesh->n_Cells();
  std::vector<double> grad_x_1, grad_x_2;
  grad_x_1.resize(n_Cell, 0.0);
  grad_x_2.resize(n_Cell, 0.0);
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    for (unsigned int j = 0; j < (it->second).size(); j++)
    {
      Face & face = (it->second)[j];
      long int node_id1 = face.node_id1();
      long int node_id2 = face.node_id2();
      long int cell_id1 = face.cell_id1();
      long int cell_id2 = face.cell_id2();

      Point point_1 = node_set[node_id1 - 1].point();
      Point point_2 = node_set[node_id2 - 1].point();
      Vec3d face_vec = point_2 - point_1;
      Vec3d face_normal;
      face_normal.x() = face_vec.y();
      face_normal.y() = -face_vec.x();
      face_normal.z() = 0.0;

      double face_x = (point_2.x() + point_1.x()) * 0.5;
      if (cell_id1 > 0)
      {
        grad_x_1[cell_id1-1] += face_x * face_normal.x() / cell_set[cell_id1-1].volume();
        grad_x_2[cell_id1-1] += face_x * face_normal.y() / cell_set[cell_id1-1].volume();
      }
      if (cell_id2 > 0)
      {
        grad_x_1[cell_id2-1] -= face_x * face_normal.x() / cell_set[cell_id2-1].volume();
        grad_x_2[cell_id2-1] -= face_x * face_normal.y() / cell_set[cell_id2-1].volume();
      }
    }
  }
  *
  for (int i = 0; i < grad_x_1.size(); i++)
    std::cout << grad_x_1[i] << std::endl;
  for (int i = 0; i < grad_x_2.size(); i++)
    std::cout << grad_x_2[i] << std::endl;*

  // solve T_cell
  PETSC_APPCTX    appCtx;
  //Vec             T_cell;
  //Vec             res;
  SNES            snes;
  Mat             J, J_MatrixFree;
  PetscInt        n_DOF = n_Cell;

  VecCreate(PETSC_COMM_SELF, &(appCtx.T_cell));
  VecSetSizes(appCtx.T_cell, PETSC_DECIDE, n_DOF);
  VecSetFromOptions(appCtx.T_cell);
  VecDuplicate(appCtx.T_cell, &(appCtx.grad_T_x));
  VecDuplicate(appCtx.T_cell, &(appCtx.grad_T_y));
  VecDuplicate(appCtx.T_cell, &(appCtx.res));

  VecCreate(PETSC_COMM_SELF, &(appCtx.T_face));
  std::cout << "p_mesh->n_Faces() = " << p_mesh->n_Faces() << std::endl;
  VecSetSizes(appCtx.T_face, PETSC_DECIDE, p_mesh->n_Faces());
  VecSetFromOptions(appCtx.T_face);
  appCtx.p_mesh = p_mesh;

  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetFunction(snes, appCtx.res, SNESFormFunction, &appCtx);

  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOF, n_DOF, 10, NULL, &J);  // 10 is the max possible neighbor cell numbers for a cell
  ISColoring   iscoloring;
  for (std::vector<FluentTriCell>::const_iterator it = cell_set.begin() ; it != cell_set.end(); ++it)
  {
    const FluentTriCell & this_cell = *it;
    PetscInt row = this_cell.id() - 1;
    PetscInt col[1];
    PetscReal array_one[1]; array_one[0] = 1.0;
    const std::vector<long int> & nbs = this_cell.getNeighborCellIDs();

    for (int j = 0; j < nbs.size(); j++)
    {
      long int nb_cell_id = nbs[j];
      if (nb_cell_id > 0)
      {
        //std::cout << "  col = " << nb_cell_id-1 << std::endl;
        const FluentTriCell & nb_cell = cell_set.at(nbs[j]-1);
        col[0] = nb_cell_id-1;
        MatSetValues(J, 1, &row, 1, col, array_one, INSERT_VALUES);
        const std::vector<long int> & nbnbs = nb_cell.getNeighborCellIDs();
        for (int k = 0; k < nbnbs.size(); k++)
        {
          long int nbnb_cell_id = nbnbs[k];
          if (nbnb_cell_id>0)
          {
            col[0] = nbnb_cell_id-1;
            //std::cout << "    col = " << nbnbs[k]-1 << std::endl;
            MatSetValues(J,             // The target matrix, Jocobian matrix J here
                   1,             // How many rows is added to the matrix, 1 here as we add one row one time
                   &row,          // Indicies of rows. Since we add one row, it is the row number
                   1,          // How many columns to be added
                   col,           // Indicies of columns
                   array_one,     // Since we want the pattern only, we give '1' to the entry
                   INSERT_VALUES);  // Inserv value
          }
        }
      }
    }
  }
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  //MatView(J, PETSC_VIEWER_STDOUT_SELF);
  //MatView(J, PETSC_VIEWER_DRAW_WORLD);

  MatCreateSNESMF(snes, &J_MatrixFree);

  MatColoring     mc;
  MatFDColoring   fdcoloring;
  MatColoringCreate(J, &mc);
  MatColoringSetType(mc, MATCOLORINGSL);
  MatColoringSetFromOptions(mc);
  MatColoringApply(mc, &iscoloring);
  MatColoringDestroy(&mc);
  MatFDColoringCreate(J, iscoloring, &fdcoloring);
  MatFDColoringSetFunction(fdcoloring, (PetscErrorCode (*)(void))(SNESFormFunction), &appCtx);
  MatFDColoringSetFromOptions(fdcoloring);
  MatFDColoringSetUp(J, iscoloring, fdcoloring);
  ISColoringDestroy(&iscoloring);
  SNESSetJacobian(snes, // snes object
                J_MatrixFree,
                J,    // Preconditioning matrix, in this case, use the same as the Jacobian matrix
                SNESComputeJacobianDefaultColor,  // Versions 3.4.3; 3.5.2
                fdcoloring);  // user object context, has to be fdcoloring here
  ISColoringDestroy(&iscoloring);

  SNESSetFromOptions(snes);

  SNESSolve(snes, NULL, appCtx.T_cell);

  std::cout << "end of SNESSolve\n";

  // error analysis
  double error = 0.0;
  PetscScalar * T_cc;
  VecGetArray(appCtx.T_cell, &T_cc);
  for (std::vector<FluentTriCell>::const_iterator it = cell_set.begin() ; it != cell_set.end(); ++it)
  {
    const FluentTriCell & this_cell = *it;
    double xx = this_cell.centroid().x();
    double yy = this_cell.centroid().y();
    double T_ana = sin(3.1415926 * xx) * sin(3.1415926 * yy);
    double T_num = T_cc[this_cell.id()-1];
    error += (T_num - T_ana) * (T_num - T_ana) * this_cell.volume();
  }
  std::cout << "L2 ERROR norm = " << sqrt(error) << std::endl;
  VecRestoreArray(appCtx.T_cell, &T_cc);

  //VecView(appCtx.T_cell, PETSC_VIEWER_STDOUT_WORLD);
  //VecView(appCtx.res, PETSC_VIEWER_STDOUT_WORLD);

  std::cout << "1\n";
  // Linear system
  GRAD gradT[n_DOF];
  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second)[j];
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          gradT[cell_id-1].addBCContribution(0.0, 0.0);
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
          gradT[cell_id1-1].addCoef(0, cell_id1, alpha_12 * face_normal.x() / v1, alpha_12 * face_normal.y() / v1);
          gradT[cell_id1-1].size++;
          gradT[cell_id1-1].addCoef(gradT[cell_id1-1].size, cell_id2, alpha_21 * face_normal.x() / v1, alpha_21 * face_normal.y() / v1);
          // Cell 2
          gradT[cell_id2-1].addCoef(0, cell_id2, -alpha_21 * face_normal.x() / v2, -alpha_21 * face_normal.y() / v2);
          gradT[cell_id2-1].size++;
          gradT[cell_id2-1].addCoef(gradT[cell_id2-1].size, cell_id1, -alpha_12 * face_normal.x() / v2, -alpha_12 * face_normal.y() / v2);
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  std::cout << "2\n";
  VecGetArray(appCtx.T_cell, &T_cc);
  std::vector<double> gradTx_test, gradTy_test;
  gradTx_test.resize(n_DOF, 0.0); gradTy_test.resize(n_DOF, 0.0);
  for (int i = 0; i < n_DOF; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      long int idx = gradT[i].cell_id[j];
      double coef_x = gradT[i].coef_x[j];
      double coef_y = gradT[i].coef_y[j];
      if (idx > 0)
      {
        double T = T_cc[idx-1];
        gradTx_test[i] += coef_x * T;
        gradTy_test[i] += coef_y * T;
      }
    }
  }

  std::cout << "3\n";
  //KSP Part
  Vec       T_KSP, b_KSP;
  Mat       M_KSP;
  KSP       ksp;
  VecDuplicate(appCtx.T_cell, &T_KSP);
  VecDuplicate(appCtx.T_cell, &b_KSP);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n_DOF, n_DOF, 10, NULL, &M_KSP);  // 10 is the max possible neighbor cell numbers for a cell
  MatSetOption(M_KSP, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, M_KSP, M_KSP);
  KSPSetFromOptions(ksp);

  for (std::map<int, std::vector<Face> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    switch (zone)
    {
      case 5:
      case 2:
      case 3:
      case 4:
      {
        for (unsigned int j = 0; j < (it->second).size(); j++)
        {
          Face & face = (it->second).at(j);
          long int cell_id1 = face.cell_id1();
          long int cell_id2 = face.cell_id2();

          long int cell_id = (cell_id1 > 0) ? cell_id1 : cell_id2;
          const FluentTriCell & cell = cell_set.at(cell_id-1);
          double distance = 0.0;

          if (zone == 5)        distance = 1.0 - cell.centroid().y();   // top
          else if (zone == 2)   distance = 1.0 - cell.centroid().x();   // right
          else if (zone == 3)   distance = cell.centroid().y();         // bottom
          else                  distance = cell.centroid().x();         // left

          PetscInt row = cell_id - 1;
          PetscInt col[1]; col[0] = row;
          double val[1];
          val[0] = face.area() / distance;
          MatSetValues(M_KSP, 1, &row, 1, col, val, ADD_VALUES);
        }
      }
      break;

      case 7:
      {
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

          col[0] = r1; val[0] = area_f / distance; nval[0] = -val[0];
          MatSetValues(M_KSP, 1, &r1, 1, col, val, ADD_VALUES);   // cell 1 diag
          MatSetValues(M_KSP, 1, &r2, 1, col, nval, ADD_VALUES);  // cell 2 off-diag
          col[0] = r2;
          MatSetValues(M_KSP, 1, &r1, 1, col, nval, ADD_VALUES);  // cell 1 off-diag
          MatSetValues(M_KSP, 1, &r2, 1, col, val, ADD_VALUES);   // cell 2 diag

          int size1 = gradT[cell_id1-1].size;
          if (size1 > 3) {std::cerr<<"size1 > 3\n"; exit(1);}
          PetscInt col1[size1+1]; double val1[size1+1], nval1[size1+1];
          for (int i = 0; i < size1+1; i++)
          {
            col1[i] = gradT[cell_id1-1].cell_id[i]-1;
            val1[i] = -gradT[cell_id1-1].coef_x[i] * n2.x() - gradT[cell_id1-1].coef_y[i] * n2.y();
            val1[i] *= (area_f * alpha_12);
            nval1[i] = -val1[i];
          }
          MatSetValues(M_KSP, 1, &r1, size1+1, col1, val1, ADD_VALUES); // cell 1 off-diag
          MatSetValues(M_KSP, 1, &r2, size1+1, col1, nval1, ADD_VALUES); // cell 2 off-diag
          int size2 = gradT[cell_id2-1].size;
          if (size2 > 3) {std::cerr<<"size2 > 3\n"; exit(1);}
          PetscInt col2[size2+1]; double val2[size2+1], nval2[size2+1];
          for (int i = 0; i < size2+1; i++)
          {
            col2[i] = gradT[cell_id2-1].cell_id[i]-1;
            val2[i] = -gradT[cell_id2-1].coef_x[i] * n2.x() - gradT[cell_id2-1].coef_y[i] * n2.y();
            val2[i] *= (area_f * alpha_21);
            nval2[i] = -val2[i];
          }
          MatSetValues(M_KSP, 1, &r1, size2+1, col2, val2, ADD_VALUES);
          MatSetValues(M_KSP, 1, &r2, size2+1, col2, nval2, ADD_VALUES);
        }
      }
      break;

      default:
        std::cerr << "ERROR" << std::endl;
    }
  }
  MatAssemblyBegin(M_KSP, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_KSP, MAT_FINAL_ASSEMBLY);

  //MatView(M_KSP, PETSC_VIEWER_STDOUT_SELF);

  std::cout << "4\n";
  PetscScalar * rhs;
  VecGetArray(b_KSP, &rhs);
  for (int i = 0; i < cell_set.size(); i++)
  {
    const FluentTriCell & cell = cell_set.at(i);
    double xx = cell.centroid().x();
    double yy = cell.centroid().y();
    double pi = 3.1415926;
    double qq = 2.0 * pi * pi * sin(pi * xx) * sin(pi * yy);
    rhs[i] = qq * cell.volume();
  }
  VecRestoreArray(b_KSP, &rhs);

  KSPSolve(ksp, b_KSP, T_KSP);

  std::cout << "5\n";
  //p_mesh->WriteVTUFile();
  FILE * ptr_File;
  ptr_File = fopen("output/FluentTwoDMesh.vtu", "w");
  p_mesh->writeMesh(ptr_File);
  std::ostringstream out_string_stream;
  out_string_stream << "      <CellData>" << "\n";

  std::cout << "5.3\n";

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
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"T_cell\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << T_cc[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";

  std::cout << "5.3\n";

  PetscScalar * TT_KSP;
  VecGetArray(T_KSP, &TT_KSP);
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"T_KSP\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << TT_KSP[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  VecRestoreArray(T_KSP, &TT_KSP);
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"gradTx\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << gradTx_test[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"gradTy\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << gradTy_test[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  // CELL DATA (temperature)
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
  {
    double x = cell_set[i].centroid().x(); double y = cell_set[i].centroid().y();
    out_string_stream << "          " << sin(2.0 * 3.14159265 * x) * cos(3.14159265 * y) << "\n";
  }
  out_string_stream << "        </DataArray>" << "\n";
  out_string_stream << "      </CellData>" << "\n";
  // POINT DATA
  out_string_stream << "      <PointData>" << "\n";
  // NODE ID
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"Node_ID\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < node_set.size(); i++)
    out_string_stream << "          " << node_set[i].id() << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  out_string_stream << "      </PointData>" << "\n";

  std::cout << "5.5\n";

  fprintf(ptr_File, "%s", out_string_stream.str().c_str());

  p_mesh->finishFile(ptr_File);
  fclose(ptr_File);
  delete p_mesh;

  VecRestoreArray(appCtx.T_cell, &T_cc);

  std::cout << "6\n";

  VecDestroy(&(appCtx.T_cell)); VecDestroy(&(appCtx.res)); VecDestroy(&(appCtx.T_face));
  VecDestroy(&(appCtx.grad_T_x)); VecDestroy(&(appCtx.grad_T_y));
  //PetscFree(appCtx.T_faces); PetscFree(appCtx.cell_vol); PetscFree(appCtx.cell_centroid); PetscFree(appCtx.grad_T_cell);
  MatDestroy(&J); MatDestroy(&J_MatrixFree);
  SNESDestroy(&snes);
  MatFDColoringDestroy(&fdcoloring);

  std::cout << "7\n";

  //KSP
  VecDestroy(&T_KSP); VecDestroy(&b_KSP); MatDestroy(&M_KSP);
  KSPDestroy(&ksp);
  PetscFinalize();

  return 0;
}
*/
