#include <cstdio>
#include <iostream>
#include <petscsnes.h>

#include "inc_problems.h"
//#include "FluentTwoDMesh.h"

std::map<int, bool>   has_U_BC, has_V_BC, has_p_BC;
std::map<int, double> U_BC, V_BC, p_BC;
std::map<int, bool>   is_interior_zone;
double DT, VISC;

void setupBoundaryConditions(FluentTwoDMesh * p_mesh)
{
  std::map<int, std::vector<Face*> > & face_zone_map = p_mesh->getFaceZoneMap();
  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;

    char type_opt[64], bc_type[64];
    PetscBool type_flg = PETSC_FALSE;
    snprintf(type_opt, sizeof(type_opt), "-bc_%d_type", zone);
    PetscOptionsGetString(NULL, NULL, type_opt, bc_type, sizeof(bc_type), &type_flg);
    if (!type_flg)
    {
      std::cout << "ERROR: mesh has zone " << zone << " but the input file does not define '"
                << type_opt << "'." << std::endl;
      exit(1);
    }
    std::string type_str(bc_type);

    if (type_str == "interior")
    {
      is_interior_zone[zone] = true;
      continue;
    }
    is_interior_zone[zone] = false;

    if (type_str == "velocity_inlet")
    {
      char u_opt[64], v_opt[64];
      snprintf(u_opt, sizeof(u_opt), "-bc_%d_u", zone);
      snprintf(v_opt, sizeof(v_opt), "-bc_%d_v", zone);
      double u_val = 0.0, v_val = 0.0;
      PetscBool u_flg = PETSC_FALSE, v_flg = PETSC_FALSE;
      PetscOptionsGetReal(NULL, NULL, u_opt, &u_val, &u_flg);
      PetscOptionsGetReal(NULL, NULL, v_opt, &v_val, &v_flg);
      if (!u_flg || !v_flg)
      {
        std::cout << "ERROR: zone " << zone << " is 'velocity_inlet' but '" << u_opt << "' and/or '" << v_opt
                  << "' are not specified." << std::endl;
        exit(1);
      }
      has_U_BC[zone] = true;  U_BC[zone] = u_val;
      has_V_BC[zone] = true;  V_BC[zone] = v_val;
      has_p_BC[zone] = false;
    }
    else if (type_str == "wall")
    {
      has_U_BC[zone] = true;  U_BC[zone] = 0.0;
      has_V_BC[zone] = true;  V_BC[zone] = 0.0;
      has_p_BC[zone] = false;
    }
    else if (type_str == "horizontal_slip_wall")
    {
      has_U_BC[zone] = false;
      has_V_BC[zone] = true;  V_BC[zone] = 0.0;
      has_p_BC[zone] = false;
    }
    else if (type_str == "vertical_slip_wall")
    {
      has_U_BC[zone] = true;  U_BC[zone] = 0.0;
      has_V_BC[zone] = false;
      has_p_BC[zone] = false;
    }
    else if (type_str == "pressure_outlet")
    {
      char p_opt[64];
      snprintf(p_opt, sizeof(p_opt), "-bc_%d_p", zone);
      double p_val = 0.0;
      PetscBool p_flg = PETSC_FALSE;
      PetscOptionsGetReal(NULL, NULL, p_opt, &p_val, &p_flg);
      if (!p_flg)
      {
        std::cout << "ERROR: zone " << zone << " is 'pressure_outlet' but '" << p_opt << "' is not specified." << std::endl;
        exit(1);
      }
      has_U_BC[zone] = false;
      has_V_BC[zone] = false;
      has_p_BC[zone] = true;  p_BC[zone] = p_val;
    }
    else
    {
      std::cout << "ERROR: zone " << zone << " has unknown -bc_" << zone << "_type '" << type_str << "'. Supported types: "
                << "interior, velocity_inlet, wall, horizontal_slip_wall, vertical_slip_wall, pressure_outlet." << std::endl;
      exit(1);
    }
  }

  /*
   * Historical alternate B.C. configuration, kept for reference while cavity support
   * is being worked out (relocated here, verbatim, from an earlier version of
   * inc_problems.h):
   *
   * static std::map<int, bool> has_U_BC = {{1, true}, {2, true}, {3, false}, {4, true}, {5, true}};
   * static std::map<int, double> U_BC = {{1, 1.0}, {2, 0.0}, {4, 0.0}, {5, 0.0}};
   *
   * has_U_BC[2] = true; has_U_BC[3] = true; has_U_BC[4] = true; has_U_BC[5] = true;
   * has_V_BC[2] = true; has_V_BC[3] = true; has_V_BC[4] = true; has_V_BC[5] = true;
   * has_p_BC[2] = false; has_p_BC[3] = false; has_p_BC[4] = false; has_p_BC[5] = false;
   * U_BC[2] = 0.0; U_BC[3] = 0.0; U_BC[4] = 0.0; U_BC[5] = 1.0;
   * V_BC[2] = 0.0; V_BC[3] = 0.0; V_BC[4] = 0.0; V_BC[5] = 0.0;
   * p_BC[2] = -1e6; p_BC[3] = -1e6; p_BC[4] = -1e6; p_BC[5] = -1e6;
   */
}

void
StreamFunctionNode::computeStreamFunction(double value)
{
  if (!_has_computed)
  {
    _stream_function_value = value;
    _has_computed = true;

    std::vector<Face*> & connected_faces = _node->getConnectedFace();
    for (unsigned int i = 0; i < connected_faces.size(); i++)
    {
      long int neighbor_node_id = connected_faces[i]->neighbor_node_id(_node->id());
      StreamFunctionNode * nb_SF_Node = _parent_vec.at(neighbor_node_id-1);
      if (!nb_SF_Node->hasComputed())
      {
        long int face_id = connected_faces[i]->id();
        double F_value = _F_face[face_id];
        double next_node_value = 0.0;

        if (_node->id() == connected_faces[i]->node_id1())
          next_node_value = _stream_function_value + F_value;
        else
          next_node_value = _stream_function_value - F_value;

        nb_SF_Node->computeStreamFunction(next_node_value);
      }
    }
  }
}

void updateAdvectionOperator(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec b_USTAR, Vec b_VSTAR, Mat M_USTAR, Mat M_VSTAR)
{
  PetscScalar * ff;
  PetscScalar * bb_ustar, * bb_vstar;
  VecGetArray(F_face_star, &ff);
  VecGetArray(b_USTAR, &bb_ustar);
  VecGetArray(b_VSTAR, &bb_vstar);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
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
        double fface = ff[face->id()];

        PetscInt col[1]; double val[1];
        PetscInt r1 = cell_id1 - 1;
        val[0] = fface; col[0] = r1;

        if (has_U_BC[zone])
          bb_ustar[cell_id1-1] -= fface * U_BC[zone];
        else
          MatSetValues(M_USTAR, 1, &r1, 1, col, val, ADD_VALUES);

        if (has_V_BC[zone])
          bb_vstar[cell_id1-1] -= fface * V_BC[zone];
        else
          MatSetValues(M_VSTAR, 1, &r1, 1, col, val, ADD_VALUES);
      }
    }
    else
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        double GAMMA = 0.5;

        Face * face = (it->second).at(j);
        long int face_id = face->id();
        long int cell_id1 = face->cell_id1();
        long int cell_id2 = face->cell_id2();

        const FluentTriCell & cell_1 = cell_set.at(cell_id1-1);
        const FluentTriCell & cell_2 = cell_set.at(cell_id2-1);
        double v1 = cell_1.volume();
        double v2 = cell_2.volume();

        double area_f = face->area();
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
  }
  VecRestoreArray(F_face_star, &ff);
  VecRestoreArray(b_USTAR, &bb_ustar);
  VecRestoreArray(b_VSTAR, &bb_vstar);

  MatAssemblyBegin(M_USTAR, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_USTAR, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(M_VSTAR, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M_VSTAR, MAT_FINAL_ASSEMBLY);
}

void updateMassVeclocities(FluentTwoDMesh * p_mesh, Vec u_STAR, Vec v_STAR, Vec F_0f_star, Vec b_p)
{
  // Update "mass velocity" flux F_0f_star and use it as pressure equation's right-hand-side
  PetscScalar * f0f;
  PetscScalar * uu_star;
  PetscScalar * vv_star;
  PetscScalar * bb_p;
  VecGetArray(F_0f_star, &f0f);
  VecGetArray(u_STAR, &uu_star);
  VecGetArray(v_STAR, &vv_star);
  VecGetArray(b_p, &bb_p);

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
        double u_face = has_U_BC[zone] ? U_BC[zone] : uu_star[cell_id1-1];
        double v_face = has_V_BC[zone] ? V_BC[zone] : vv_star[cell_id1-1];
        double fface = u_face * face_normal.x() + v_face * face_normal.y();

        f0f[face->id()] = fface;
        bb_p[cell_id1-1] += fface / DT;
      }
    }
    else
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        Face * face = (it->second).at(j);
        double r = face->distance_ratio();

        long int cell_id1 = face->cell_id1();
        long int cell_id2 = face->cell_id2();

        double u_face = uu_star[cell_id1-1] * (1.0 - r) + uu_star[cell_id2-1] * r;
        double v_face = vv_star[cell_id1-1] * (1.0 - r) + vv_star[cell_id2-1] * r;

        Vec3d face_normal = face->faceNormal();

        double fface = u_face * face_normal.x() + v_face * face_normal.y();
        f0f[face->id()] = fface;

        bb_p[cell_id1-1] += fface / DT;
        bb_p[cell_id2-1] -= fface / DT;
      }
    }

    // anchor point
    //bb_p[0] = 0.0;
  }
  VecRestoreArray(F_0f_star, &f0f);
  VecRestoreArray(u_STAR, &uu_star);
  VecRestoreArray(v_STAR, &vv_star);
  VecRestoreArray(b_p, &bb_p);
}

void updateFfaceStar(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec F_0f_star, Vec p, Vec gradP_x, Vec gradP_y)
{
  PetscScalar * ff, * f0f, * pp;
  PetscScalar * grad_p_x, * grad_p_y;
  VecGetArray(F_face_star, &ff);  VecGetArray(F_0f_star, &f0f); VecGetArray(p, &pp);
  VecGetArray(gradP_x, &grad_p_x);  VecGetArray(gradP_y, &grad_p_y);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  std::map<int, std::vector<Face*> > & face_zone_map = p_mesh->getFaceZoneMap();
  for (std::map<int, std::vector<Face*> >::iterator it = face_zone_map.begin(); it != face_zone_map.end(); ++it)
  {
    int zone = it->first;
    if (!is_interior_zone[zone])
    {
      for (unsigned int j = 0; j < (it->second).size(); j++)
      {
        // Copy from F_0f_star (as boundary conditions have been setup there)
        Face * face = (it->second).at(j);
        ff[face->id()] = f0f[face->id()];
      }
    }
    else
    {
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

        double cross_diff = alpha_12 * (grad_p_x[cell_id1-1] * n2.x() + grad_p_y[cell_id1-1] * n2.y())
                          + alpha_21 * (grad_p_x[cell_id2-1] * n2.x() + grad_p_y[cell_id2-1] * n2.y());
        double pressure_correction = area_f * ((pp[cell_id2-1] - pp[cell_id1-1]) / distance + cross_diff);

        ff[face->id()] = f0f[face->id()] - DT * pressure_correction;
      }
    }
  }
  VecRestoreArray(F_face_star, &ff); VecRestoreArray(F_0f_star, &f0f); VecRestoreArray(p, &pp);
  VecRestoreArray(gradP_x, &grad_p_x); VecRestoreArray(gradP_y, &grad_p_y);
}

void evaluatePressureGradientValues(FluentTwoDMesh * p_mesh, Vec p, Vec gradP_x, Vec gradP_y, GRAD * grad_p)
{
  VecSet(gradP_x, 0.0); VecSet(gradP_y, 0.0);
  PetscScalar * pp, * grad_p_x, * grad_p_y;
  VecGetArray(p, &pp); VecGetArray(gradP_x, &grad_p_x); VecGetArray(gradP_y, &grad_p_y);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();

  for(std::vector<FluentTriCell>::iterator it = cell_set.begin(); it != cell_set.end(); ++it)
  {
    long int cell_id = it->id();
    GRAD & grad_p_cell = grad_p[cell_id-1];
    for (int i = 0; i < grad_p_cell.size+1; i++)
    {
      long int cell_id_i = grad_p_cell.cell_id[i];
      double pi = pp[cell_id_i-1];
      grad_p_x[cell_id-1] += grad_p_cell.coef_x[i] * pi;
      grad_p_y[cell_id-1] += grad_p_cell.coef_y[i] * pi;
    }
    grad_p_x[cell_id-1] += grad_p_cell.bc_x;
    grad_p_y[cell_id-1] += grad_p_cell.bc_y;
  }
  VecRestoreArray(gradP_x, &grad_p_x); VecRestoreArray(gradP_y, &grad_p_y); VecRestoreArray(p, &pp);
}

void updatePressureGradientAsSource(FluentTwoDMesh * p_mesh, Vec b_USTAR, Vec b_VSTAR, Vec gradP_x, Vec gradP_y)
{
  PetscScalar * bb_ustar, * bb_vstar, * grad_p_x, * grad_p_y;
  VecGetArray(b_USTAR, &bb_ustar); VecGetArray(b_VSTAR, &bb_vstar);
  VecGetArray(gradP_x, &grad_p_x); VecGetArray(gradP_y, &grad_p_y);

  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();

  for(std::vector<FluentTriCell>::iterator it = cell_set.begin(); it != cell_set.end(); ++it)
  {
    long int cell_id = it->id();
    bb_ustar[cell_id-1] -= grad_p_x[cell_id-1] * it->volume();
    bb_vstar[cell_id-1] -= grad_p_y[cell_id-1] * it->volume();
  }

  VecRestoreArray(b_USTAR, &bb_ustar);  VecRestoreArray(b_VSTAR, &bb_vstar);
  VecRestoreArray(gradP_x, &grad_p_x);  VecRestoreArray(gradP_y, &grad_p_y);
}


/**************
 * output
 **************/
void writeOutputFile(int time_step, FluentTwoDMesh * p_mesh, Vec u, Vec v, Vec p, Vec F_face_star, bool computeStreamFunction)
{
  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();
  std::vector<Node*> & node_set = p_mesh->getNodeSet();

  std::string file_name = "output/output_" + std::to_string(time_step) + ".vtu";
  FILE * ptr_File;
  ptr_File = fopen(file_name.c_str(), "w");
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

  PetscScalar * uu;
  VecGetArray(u, &uu);
  out_string_stream << "        <DataArray type=\"Float32\" Name=\"u_star\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < cell_set.size(); i++)
    out_string_stream << "          " << uu[i] << "\n";
  out_string_stream << "        </DataArray>" << "\n";
  VecRestoreArray(u, &uu);

  PetscScalar * vv;
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
    out_string_stream << "          " << node_set[i]->id() << "\n";
  out_string_stream << "        </DataArray>" << "\n";

  // Stream function
  // Assuming it is expensive to evaluate streamline function, just compute the last time step
  /**** Stream Function ****/
  if (computeStreamFunction)
  {
    std::vector<StreamFunctionNode *> SFNodeVec(node_set.size(), NULL);

    PetscScalar * ff;
    VecGetArray(F_face_star, &ff);
    for (int i = 0; i < node_set.size(); i++)
    {
      StreamFunctionNode * SFNode = new StreamFunctionNode(node_set.at(i), SFNodeVec, ff);
      SFNodeVec[i] = SFNode;
    }

    SFNodeVec.at(0)->computeStreamFunction(0.0);

    VecRestoreArray(F_face_star, &ff);

    out_string_stream << "        <DataArray type=\"Float32\" Name=\"StreamFunction\" format=\"ascii\">" << "\n";
    for(unsigned int i = 0; i < node_set.size(); i++)
      out_string_stream << "          " << SFNodeVec[i]->getSFValue() << "\n";
    out_string_stream << "        </DataArray>" << "\n";
    for (int i = 0; i < SFNodeVec.size(); i++)
    {
      delete SFNodeVec.at(i);
    }
  }

  out_string_stream << "      </PointData>" << "\n";

  fprintf(ptr_File, "%s", out_string_stream.str().c_str());

  p_mesh->finishFile(ptr_File);
  fclose(ptr_File);
}
