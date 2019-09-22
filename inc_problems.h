#ifndef INC_PROBLEMS_H
#define INC_PROBLEMS_H

#include "FluentTwoDMesh.h"

static double DT = 0.01;
static double VISC = 1.0;

struct GRAD;
void updateAdvectionOperator(FluentTwoDMesh * p_mesh, Vec F_face_star, Mat M_USTAR, Mat M_VSTAR);
void updateMassVeclocities(FluentTwoDMesh * p_mesh, Vec u_STAR, Vec v_STAR, Vec F_0f_star, Vec b_p);
void updateFfaceStar(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec F_0f_star, Vec p, GRAD * grad_u_star);
void updatePressureGradientAsSource(FluentTwoDMesh * p_mesh, Vec p, Vec p_src_x, Vec p_src_y);

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

#endif // INC_PROBLEMS_H
