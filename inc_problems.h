#ifndef INC_PROBLEMS_H
#define INC_PROBLEMS_H

#include "FluentTwoDMesh.h"
#include <map>

/* Setup B.C. */
static std::map<int, bool> has_U_BC;
static std::map<int, bool> has_V_BC;
static std::map<int, bool> has_p_BC;
static std::map<int, double> U_BC;
static std::map<int, double> V_BC;
static std::map<int, double> p_BC;
/*
has_U_BC[2] = true; has_U_BC[3] = true; has_U_BC[4] = true; has_U_BC[5] = true;
has_V_BC[2] = true; has_V_BC[3] = true; has_V_BC[4] = true; has_V_BC[5] = true;
has_p_BC[2] = false; has_p_BC[3] = false; has_p_BC[4] = false; has_p_BC[5] = false;
U_BC[2] = 0.0; U_BC[3] = 0.0; U_BC[4] = 0.0; U_BC[5] = 1.0;
V_BC[2] = 0.0; V_BC[3] = 0.0; V_BC[4] = 0.0; V_BC[5] = 0.0;
p_BC[2] = -1e6; p_BC[3] = -1e6; p_BC[4] = -1e6; p_BC[5] = -1e6;*/
/* End of B.C. setting */

static double DT = 0.05;
static double VISC = 0.01;

class StreamFunctionNode
{
public:
  StreamFunctionNode(Node * node, std::vector<StreamFunctionNode *> & parent_vec, double * F_face) :
  _node(node), _has_computed(false), _stream_function_value(0.0), _parent_vec(parent_vec), _F_face(F_face)
  {}
  ~StreamFunctionNode()
  {}

  void computeStreamFunction(double value);
  double getSFValue() { return _stream_function_value; }
  void setSFValue(double val) { _stream_function_value = val; }
  bool hasComputed() { return _has_computed; }

protected:
  Node * _node;
  bool _has_computed;
  double _stream_function_value;
  std::vector<StreamFunctionNode *> & _parent_vec;
  double * _F_face;
};

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
