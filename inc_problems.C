#include <iostream>
#include <petscsnes.h>

#include "inc_problems.h"
#include "FluentTwoDMesh.h"

void updateAdvectionOperator(FluentTwoDMesh * p_mesh, Vec F_face_star, Mat M_USTAR, Mat M_VSTAR)
{
  PetscScalar * ff;
  VecGetArray(F_face_star, &ff);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  std::map<int, std::vector<Face> > & face_zone_map = p_mesh->getFaceZoneMap();
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
}

void updateMassVeclocities(FluentTwoDMesh * p_mesh, Vec u_STAR, Vec v_STAR, Vec F_0f_star, Vec b_p)
{
  // Update "mass velocity" flux F_0f_star
  PetscScalar * f0f;
  PetscScalar * uu_star;
  PetscScalar * vv_star;
  PetscScalar * bb_p;
  VecGetArray(F_0f_star, &f0f);
  VecGetArray(u_STAR, &uu_star);
  VecGetArray(v_STAR, &vv_star);
  VecGetArray(b_p, &bb_p);

  std::map<int, std::vector<Face> > & face_zone_map = p_mesh->getFaceZoneMap();
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
}

void updateFfaceStar(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec F_0f_star, Vec p, GRAD * grad_u_star)
{
  PetscScalar * ff;
  PetscScalar * pp;
  PetscScalar * f0f;
  VecGetArray(F_face_star, &ff);
  VecGetArray(F_0f_star, &f0f);
  VecGetArray(p, &pp);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  std::map<int, std::vector<Face> > & face_zone_map = p_mesh->getFaceZoneMap();
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
  //std::cout << "End loop: before restore vec" << std::endl;
  VecRestoreArray(F_face_star, &ff);
  VecRestoreArray(F_0f_star, &f0f);
  VecRestoreArray(p, &pp);
  //std::cout << "End loop" << std::endl;
}