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

// Optional per-zone spatially-varying replacement for a constant U_BC[zone], selected
// via "-bc_<zone>_u_profile parabolic" (see setupBoundaryConditions() below). Only u
// supports a profile -- v remains a plain constant, as today.
struct ParabolicProfile { double y_min, y_max, u_max; };
extern std::map<int, bool> has_U_profile;
extern std::map<int, ParabolicProfile> U_profile;

// Returns the Dirichlet u-value to use for `face` in `zone` -- either the constant
// U_BC[zone], or, when has_U_profile[zone], the parabolic profile evaluated at the
// face's midpoint y (u_max * (1 - eta^2), eta = normalized position in [y_min, y_max]).
// Callers must already know has_U_BC[zone] is true.
double evaluateFaceU_BC(int zone, Face * face, FluentTwoDMesh * p_mesh);

// Temperature B.C. map: optional per boundary zone, via "-bc_<zone>_T <value>" (see
// setupBoundaryConditions() below). A zone with no such option is adiabatic
// (zero-gradient), the same default a wall's u/v gets when omitted.
extern std::map<int, bool>   has_T_BC;
extern std::map<int, double> T_BC;
// Thermal diffusivity ("-diff_t"), Gamma_T in DEB08's Eq. 2. Whether this should be
// 1/(Re*Pr) (forced convection) or 1/sqrt(Ra*Pr) (natural convection) is up to the
// caller, mirroring how VISC is handed -visc directly rather than derived from Re here.
extern double DIFF_T;
// Set via "-buoyancy true/false" (default false). Gates whether updateBuoyancySource()
// is called at all -- T is always solved regardless (see cxx_main.C), so with this off
// T is a pure passive scalar with no feedback into momentum.
extern bool BUOYANCY_ON;

// Set via "-newton true/false" (default false). Selects the monolithic Newton-Krylov
// solver (newton_solver.C) in place of the default segregated Picard/SIMPLE time loop
// in cxx_main.C. The two solve the same discrete equations -- see
// docs/newton_krylov_formulation.tex -- so with -buoyancy false they must agree to
// solver tolerance; that equivalence is the Newton path's main regression test.
extern bool NEWTON_ON;

// Whether a given mesh zone is an interior zone (shared by two cells) rather than a
// boundary zone (one cell, possibly a Dirichlet B.C.) -- populated by
// setupBoundaryConditions() below from each zone's "-bc_<zone>_type", so callers no
// longer need to hardcode a specific zone number (e.g. "7") as "the" interior zone.
extern std::map<int, bool> is_interior_zone;

// True when no zone provides a pressure Dirichlet B.C. (e.g. the lid-driven cavity,
// which is walls all the way around). The pressure equation is otherwise only defined
// up to an additive constant, so main() and updateMassVeclocities() anchor cell 0's
// pressure to zero when this is set -- see the matching MatZeroRows()/bb_p[0]=0.0.
extern bool PRESSURE_ANCHOR_NEEDED;

// Reads "-bc_<zone>_type <type>" for every zone in p_mesh (plus per-type value
// options: "-bc_<zone>_u"/"-bc_<zone>_v" for velocity_inlet, "-bc_<zone>_p" for
// pressure_outlet) and populates the B.C. maps and is_interior_zone above. Supported
// types:
//   interior              - shared by two cells, no B.C. (see is_interior_zone above)
//   velocity_inlet        - u, v both Dirichlet. -bc_<zone>_v is REQUIRED. -bc_<zone>_u
//                            is REQUIRED unless "-bc_<zone>_u_profile parabolic" is given
//                            instead, in which case -bc_<zone>_u_ymin/_ymax/_umax
//                            (all REQUIRED) instead give a parabolic u(y) profile --
//                            see evaluateFaceU_BC() above
//   wall                  - u, v Dirichlet; -bc_<zone>_u/_v are OPTIONAL and default to
//                            0 (no-slip), so a moving wall (e.g. a cavity's lid) is just
//                            a wall with an explicit non-zero -bc_<zone>_u/_v
//   horizontal_slip_wall  - v = 0 only (u left free); e.g. a horizontal tunnel wall
//   vertical_slip_wall    - u = 0 only (v left free); e.g. a vertical tunnel wall
//   pressure_outlet       - p Dirichlet, value from -bc_<zone>_p
// A generic (non-axis-aligned) slip wall isn't supported yet -- only the two special
// cases above. Exits with a clear error if a mesh zone has no matching "-bc_<zone>_type"
// or the type name/required values are missing.
//
// Independently of the above, every non-interior zone may also specify "-bc_<zone>_T
// <value>" for a Dirichlet temperature B.C.; if omitted, the zone is adiabatic.
//
// Also sets PRESSURE_ANCHOR_NEEDED (see above) once every zone has been processed.
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
// Same structure as updateAdvectionOperator() above, but for the single scalar field T.
void updateEnergyAdvectionOperator(FluentTwoDMesh * p_mesh, Vec F_face, Vec b_T, Mat M_T);
// Adds the Boussinesq buoyancy source (S_v = T, DEB08 section 4's Ra/Pr-nondimensionalized
// form) into b_VSTAR. Callers gate this on BUOYANCY_ON themselves, so the toggle stays
// visible at the call site in cxx_main.C rather than hidden inside this function.
void updateBuoyancySource(FluentTwoDMesh * p_mesh, Vec T, Vec b_VSTAR);
void updateMassVeclocities(FluentTwoDMesh * p_mesh, Vec u_STAR, Vec v_STAR, Vec F_0f_star, Vec b_p);
void updateFfaceStar(FluentTwoDMesh * p_mesh, Vec F_face_star, Vec F_0f_star, Vec p, Vec gradP_x, Vec gradP_y);
void evaluatePressureGradientValues(FluentTwoDMesh * p_mesh, Vec p, Vec gradP_x, Vec gradP_y, GRAD * grad_p);
//void updatePressureGradientAsSource(FluentTwoDMesh * p_mesh, Vec p, Vec p_src_x, Vec p_src_y, GRAD * grad_p);
void updatePressureGradientAsSource(FluentTwoDMesh * p_mesh, Vec b_USTAR, Vec b_VSTAR, Vec grad_p_x, Vec grad_p_y);

// Writes "<output_prefix>_<time_step>.vtu", e.g. output_prefix "output/cylinder_out"
// produces "output/cylinder_out_0.vtu", "output/cylinder_out_1.vtu", etc.
void writeOutputFile(const std::string & output_prefix, int time_step, FluentTwoDMesh * p_mesh, Vec u, Vec v, Vec p, Vec F_face_star, bool computeStreamFunction);

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
