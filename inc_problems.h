#ifndef INC_PROBLEMS_H
#define INC_PROBLEMS_H

#include "FluentTwoDMesh.h"
#include <map>
#include <string>

/*
 * B.C. maps, DT, and VISC are runtime-configured by setupBoundaryConditions() (see
 * inc_problems.C), driven by "-bc_<zone>_type" (and per-type value) options read from
 * the -i input file -- see setupBoundaryConditions()'s declaration below for the
 * supported types. These must be `extern` (defined exactly once, in inc_problems.C):
 * both this header and inc_problems.C are compiled into separate translation units,
 * and a `static` (internal-linkage) definition here would give each one its own
 * independent copy, silently out of sync the moment they're set at runtime instead of
 * at compile time.
 */
extern std::map<int, bool>   has_U_BC, has_V_BC, has_p_BC;
extern std::map<int, double> U_BC, V_BC, p_BC;
extern double DT, VISC;

// Whether a given mesh zone is an interior zone (shared by two cells) rather than a
// boundary zone (one cell, possibly a Dirichlet B.C.) -- populated by
// setupBoundaryConditions() below from each zone's "-bc_<zone>_type", so callers no
// longer need to hardcode a specific zone number (e.g. "7") as "the" interior zone.
extern std::map<int, bool> is_interior_zone;

// Reads "-bc_<zone>_type <type>" for every zone in p_mesh (plus per-type value
// options: "-bc_<zone>_u"/"-bc_<zone>_v" for velocity_inlet, "-bc_<zone>_p" for
// pressure_outlet) and populates the B.C. maps and is_interior_zone above. Supported
// types:
//   interior              - shared by two cells, no B.C. (see is_interior_zone above)
//   velocity_inlet        - u, v both Dirichlet, values from -bc_<zone>_u/_v
//   wall                  - no-slip: u = v = 0
//   horizontal_slip_wall  - v = 0 only (u left free); e.g. a horizontal tunnel wall
//   vertical_slip_wall    - u = 0 only (v left free); e.g. a vertical tunnel wall
//   pressure_outlet       - p Dirichlet, value from -bc_<zone>_p
// A generic (non-axis-aligned) slip wall isn't supported yet -- only the two special
// cases above. Exits with a clear error if a mesh zone has no matching "-bc_<zone>_type"
// or the type name/required values are missing.
//
// This does NOT yet handle problems like the lid-driven cavity, which have no
// Dirichlet pressure B.C. at all (a pressure anchor *point* instead of a zone) --
// planned as separate follow-up work.
void setupBoundaryConditions(FluentTwoDMesh * p_mesh);

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
void updateAdvectionOperator(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec b_USTAR, Vec b_VSTAR, Mat M_USTAR, Mat M_VSTAR);
void updateMassVeclocities(FluentTwoDMesh * p_mesh, Vec u_STAR, Vec v_STAR, Vec F_0f_star, Vec b_p);
void updateFfaceStar(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec F_0f_star, Vec p, Vec gradP_x, Vec gradP_y);
void evaluatePressureGradientValues(FluentTwoDMesh * p_mesh, Vec p, Vec gradP_x, Vec gradP_y, GRAD * grad_p);
//void updatePressureGradientAsSource(FluentTwoDMesh * p_mesh, Vec p, Vec p_src_x, Vec p_src_y, GRAD * grad_p);
void updatePressureGradientAsSource(FluentTwoDMesh * p_mesh, Vec b_USTAR, Vec b_VSTAR, Vec grad_p_x, Vec grad_p_y);

void writeOutputFile(int time_step, FluentTwoDMesh * p_mesh, Vec u, Vec v, Vec p, Vec F_face_star, bool computeStreamFunction);

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
