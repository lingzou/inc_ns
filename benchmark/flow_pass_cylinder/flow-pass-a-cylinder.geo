// Flow past a cylinder, enlarged domain: 50 x 25 rectangular channel with a
// circular hole, cylinder centered in the middle of the field.
//
//   D(0,25) ------------------------------------------- C(50,25)
//      |                    top wall                       |
//   inlet                 O(25,12.5)                     outlet
//   (left)                 dia=1.0                        (right)
//      |                   bottom wall                     |
//   A(0,0) -------------------------------------------- B(50,0)
//
// Boundaries (zone ids match input/cylinder.i, cylinder_newton_Schur.i and
// heated_cylinder.i's -bc_<zone> convention -- all three read a
// flow-past-a-cylinder.msh, generated from this file via gmsh_to_fluent.py):
//   AD : inlet             -- zone 1
//   AB : bottom wall       -- zone 2
//   BC : outlet            -- zone 3
//   CD : top wall          -- zone 4
//   circle : cylinder wall -- zone 5
//
// Every resolution knob below is copied verbatim from mesh/cylinder.geo -- only
// the domain extent and the cylinder's position within it change.

Lx = 50.0;
Ly = 25.0;
Ox = 25.0;  // cylinder center, in the middle of the field
Oy = 12.5;
r  = 0.5;   // cylinder radius (diameter 1.0)

Point(1) = {0,  0,  0}; // A
Point(2) = {Lx, 0,  0}; // B
Point(3) = {Lx, Ly, 0}; // C
Point(4) = {0,  Ly, 0}; // D

Line(1) = {1, 2}; // AB - bottom wall
Line(2) = {2, 3}; // BC - outlet
Line(3) = {3, 4}; // CD - top wall
Line(4) = {4, 1}; // DA - inlet

// Cylinder: center point + 4 points on the circle, 90 degrees apart, so the
// circle is built from 4 exact circular arcs (Transfinite subdivision below
// keeps every node exactly on the circle).
Point(5) = {Ox,     Oy,     0}; // center O
Point(6) = {Ox + r, Oy,     0}; //   0 deg
Point(7) = {Ox,     Oy + r, 0}; //  90 deg
Point(8) = {Ox - r, Oy,     0}; // 180 deg
Point(9) = {Ox,     Oy - r, 0}; // 270 deg

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

Physical Curve("inlet", 1)      = {4};
Physical Curve("bottom_wall", 2) = {1};
Physical Curve("outlet", 3)     = {2};
Physical Curve("top_wall", 4)   = {3};
Physical Curve("cylinder", 5)   = {5, 6, 7, 8};

Physical Surface("domain", 6) = {1};

// 80 elements around the circle (20 per quarter-arc, uniform since the
// default Transfinite distribution on a Circle divides equal arc length).
Transfinite Curve {5, 6, 7, 8} = 21;

// Mesh size grows away from the cylinder wall at a fixed geometric ratio.
// For a sequence of element sizes h_i = h0 * ratio^i, the cumulative
// distance after i layers is d_i = h0*(ratio^i - 1)/(ratio-1), which
// rearranges to the *exact* closed form h(d) = h0 + (ratio-1)*d -- so a
// plain linear MathEval field reproduces true geometric growth at the
// layer boundaries, no per-layer boundary-layer mesh required.
h0    = Pi * (2*r) / 80; // element size at the cylinder wall (matches the 80 arc divisions)
ratio = 1.05;
growth = ratio - 1;

lc_general = 0.7; // far-field cap away from the wake (2x mesh/cylinder.geo's 0.35,
                   // coarser here since the enlarged domain pushes the inlet/outlet
                   // and top/bottom walls much further from the cylinder)
lc_wake    = 0.12; // cap inside the downstream wake corridor

Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 200;

// Radial growth away from the cylinder, capped at the general far-field size.
Field[2] = MathEval;
Field[2].F = Sprintf("Min(%g + (%g)*F1, %g)", h0, growth, lc_general);

// Wake corridor: a box downstream of the cylinder (starting at its upstream
// edge so it blends with the radial field near the cylinder), same absolute
// extent as mesh/cylinder.geo's (14.5D downstream of the cylinder center) so
// the near-wake resolution matches the base case; it now trails well short of
// the outlet since the enlarged domain adds extra buffer downstream.
Field[3] = Box;
Field[3].XMin = Ox - r;
Field[3].XMax = Ox + 14.5;
Field[3].YMin = Oy - 1.5;
Field[3].YMax = Oy + 1.5;
Field[3].Thickness = 2.0;
Field[3].VIn = lc_wake;
Field[3].VOut = lc_general;

Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

// Pure triangle mesh (no recombination), matching gmsh_to_fluent.py's
// triangle-only conversion.
Mesh.Algorithm = 6; // Frontal-Delaunay
