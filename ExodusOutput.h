#ifndef EXODUS_OUTPUT_H
#define EXODUS_OUTPUT_H

#include <petscdmplex.h>
#include "FluentTwoDMesh.h"

/*
 * Output-only DMPlex mirror of the solver's FluentTwoDMesh, used exclusively to drive
 * PETSc's ExodusII viewer. It is never touched by the solver's own assembly/KSP code.
 *
 * PETSc's DMPlex ExodusII writer requires one single-field DM per output variable
 * (VecView keys off PetscObjectGetName(vec) matching a zonal/nodal variable name that
 * was pre-registered on the viewer), rather than a single combined multi-field vector.
 */
struct ExodusOutputContext
{
  DM baseDM;                   // topology + "Cell Sets" label; written once via DMView
  DM dmU, dmV, dmP, dmT, dmSF; // single-field clones of baseDM, one per output variable
  PetscViewer viewer;
  PetscInt cStart, cEnd;  // cell point range (height 0)
  PetscInt vStart, vEnd;  // vertex point range (depth 0)
};

// Builds the DMPlex mirror from p_mesh's node coordinates and cell connectivity,
// opens filename as an ExodusII file, and writes the mesh topology once. Call once,
// right after the mesh is loaded from file.
ExodusOutputContext * createExodusOutput(FluentTwoDMesh * p_mesh, const char * filename);

// Writes one time step's vel_x/vel_y/pressure/T element (zonal) variables, and, only
// when writeStreamFunction is true, the StreamFunction nodal variable computed from
// F_face_star (mirrors the same computation writeOutputFile does for VTU).
void writeExodusStep(ExodusOutputContext * ctx, FluentTwoDMesh * p_mesh, int step, double time,
                      Vec u, Vec v, Vec p, Vec T, Vec F_face_star, bool writeStreamFunction);

// Destroys the viewer and all DMs, and frees ctx.
void closeExodusOutput(ExodusOutputContext * ctx);

#endif // EXODUS_OUTPUT_H
