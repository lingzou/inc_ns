#include <vector>
#include <exodusII.h>

#include "ExodusOutput.h"
#include "inc_problems.h"

namespace {

DM makeFieldDM(DM baseDM, const char * fieldName, PetscInt pStart, PetscInt pEnd)
{
  DM fieldDM;
  DMClone(baseDM, &fieldDM);

  PetscSection section;
  PetscSectionCreate(PETSC_COMM_WORLD, &section);
  PetscSectionSetNumFields(section, 1);
  PetscSectionSetFieldName(section, 0, fieldName);
  PetscSectionSetChart(section, pStart, pEnd);
  for (PetscInt pt = pStart; pt < pEnd; pt++)
  {
    PetscSectionSetDof(section, pt, 1);
    PetscSectionSetFieldDof(section, pt, 0, 1);
  }
  PetscSectionSetUp(section);
  DMSetLocalSection(fieldDM, section);
  PetscSectionDestroy(&section);

  return fieldDM;
}

// Writes one variable's values for one time step: builds a fresh global Vec on
// fieldDM (already sectioned/named to match a variable pre-registered on the
// viewer), copies srcValues in point order, then hands it to VecView.
void writeFieldValues(DM fieldDM, PetscInt pStart, PetscInt pEnd, const PetscScalar * srcValues,
                       const char * name, PetscViewer viewer, int step, double time)
{
  PetscSection section;
  DMGetGlobalSection(fieldDM, &section);

  Vec v;
  DMCreateGlobalVector(fieldDM, &v);
  PetscObjectSetName((PetscObject)v, name);

  PetscScalar * arr;
  VecGetArray(v, &arr);
  for (PetscInt pt = pStart; pt < pEnd; pt++)
  {
    PetscInt off;
    PetscSectionGetOffset(section, pt, &off);
    arr[off] = srcValues[pt - pStart];
  }
  VecRestoreArray(v, &arr);

  DMSetOutputSequenceNumber(fieldDM, step, (PetscReal)time);
  VecView(v, viewer);
  VecDestroy(&v);
}

} // namespace

ExodusOutputContext * createExodusOutput(FluentTwoDMesh * p_mesh, const char * filename)
{
  std::vector<Node*> & node_set = p_mesh->getNodeSet();
  std::vector<FluentTriCell> & cell_set = p_mesh->getCellSet();

  PetscInt numVertices = (PetscInt)node_set.size();
  PetscInt numCells = (PetscInt)cell_set.size();
  PetscInt numCorners = 3, spaceDim = 2;

  // Flatten into the 0-based, id-indexed arrays DMPlexCreateFromCellListPetsc expects,
  // writing by (id - 1) rather than iteration order to match the 1-based ids used
  // throughout the rest of the codebase (e.g. cell_set.at(cell_id-1)).
  std::vector<PetscReal> coords(numVertices * spaceDim);
  for (size_t i = 0; i < node_set.size(); i++)
  {
    long int idx = node_set[i]->id() - 1;
    coords[idx * spaceDim + 0] = node_set[i]->x();
    coords[idx * spaceDim + 1] = node_set[i]->y();
  }

  std::vector<PetscInt> cellConn(numCells * numCorners);
  for (size_t i = 0; i < cell_set.size(); i++)
  {
    long int idx = cell_set[i].id() - 1;
    std::vector<long int> nodeIDs = cell_set[i].getNodeIDs();
    for (int k = 0; k < numCorners; k++)
      cellConn[idx * numCorners + k] = (PetscInt)(nodeIDs[k] - 1);
  }

  ExodusOutputContext * ctx = new ExodusOutputContext();

  DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, spaceDim, numCells, numVertices, numCorners,
                                 PETSC_TRUE, cellConn.data(), spaceDim, coords.data(), &ctx->baseDM);
  PetscObjectSetName((PetscObject)ctx->baseDM, "inc_ns_mesh");

  // Required for DMView to the ExodusII viewer -- without a "Cell Sets" label (mapping
  // cells to Exodus element blocks), it segfaults inside PETSc's own writer.
  DMCreateLabel(ctx->baseDM, "Cell Sets");
  DMLabel cellSetLabel;
  DMGetLabel(ctx->baseDM, "Cell Sets", &cellSetLabel);

  DMPlexGetHeightStratum(ctx->baseDM, 0, &ctx->cStart, &ctx->cEnd);
  DMPlexGetDepthStratum(ctx->baseDM, 0, &ctx->vStart, &ctx->vEnd);
  for (PetscInt c = ctx->cStart; c < ctx->cEnd; c++)
    DMLabelSetValue(cellSetLabel, c, 1);

  ctx->dmU  = makeFieldDM(ctx->baseDM, "u_star", ctx->cStart, ctx->cEnd);
  ctx->dmV  = makeFieldDM(ctx->baseDM, "v_star", ctx->cStart, ctx->cEnd);
  ctx->dmP  = makeFieldDM(ctx->baseDM, "pressure", ctx->cStart, ctx->cEnd);
  ctx->dmSF = makeFieldDM(ctx->baseDM, "StreamFunction", ctx->vStart, ctx->vEnd);

  PetscViewerExodusIIOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &ctx->viewer);
  // Without an explicit order, PETSc's writer leaves its internal degree unset and
  // segfaults dereferencing an uninitialized node-ordering table (verified against
  // PETSc 3.25.5 source, src/dm/impls/plex/exodusii/plexexodusii2.c).
  PetscViewerExodusIISetOrder(ctx->viewer, 1);
  DMView(ctx->baseDM, ctx->viewer);

  // The Exodus/netCDF schema needs its variable table declared upfront -- VecView
  // later looks up the Vec's own name against this table, it does not infer variables
  // from a DM's section field names.
  const char * zonalNames[3] = {"u_star", "v_star", "pressure"};
  const char * nodalNames[1] = {"StreamFunction"};
  PetscViewerExodusIISetZonalVariable(ctx->viewer, 3);
  PetscViewerExodusIISetZonalVariableNames(ctx->viewer, zonalNames);
  PetscViewerExodusIISetNodalVariable(ctx->viewer, 1);
  PetscViewerExodusIISetNodalVariableNames(ctx->viewer, nodalNames);

  return ctx;
}

void writeExodusStep(ExodusOutputContext * ctx, FluentTwoDMesh * p_mesh, int step, double time,
                      Vec u, Vec v, Vec p, Vec F_face_star, bool writeStreamFunction)
{
  const PetscScalar * uu;
  const PetscScalar * vv;
  const PetscScalar * pp;
  VecGetArrayRead(u, &uu);
  VecGetArrayRead(v, &vv);
  VecGetArrayRead(p, &pp);

  writeFieldValues(ctx->dmU, ctx->cStart, ctx->cEnd, uu, "u_star", ctx->viewer, step, time);
  writeFieldValues(ctx->dmV, ctx->cStart, ctx->cEnd, vv, "v_star", ctx->viewer, step, time);
  writeFieldValues(ctx->dmP, ctx->cStart, ctx->cEnd, pp, "pressure", ctx->viewer, step, time);

  VecRestoreArrayRead(u, &uu);
  VecRestoreArrayRead(v, &vv);
  VecRestoreArrayRead(p, &pp);

  if (writeStreamFunction)
  {
    // Mirrors the stream-function computation in writeOutputFile() (inc_problems.C),
    // duplicated rather than shared so the already-validated VTU writer stays
    // completely untouched by this addition.
    std::vector<Node*> & node_set = p_mesh->getNodeSet();
    std::vector<StreamFunctionNode*> SFNodeVec(node_set.size(), NULL);

    PetscScalar * ff;
    VecGetArray(F_face_star, &ff);
    for (size_t i = 0; i < node_set.size(); i++)
      SFNodeVec[i] = new StreamFunctionNode(node_set[i], SFNodeVec, ff);
    SFNodeVec.at(0)->computeStreamFunction(0.0);
    VecRestoreArray(F_face_star, &ff);

    std::vector<PetscScalar> sfValues(node_set.size());
    for (size_t i = 0; i < node_set.size(); i++)
    {
      sfValues[i] = SFNodeVec[i]->getSFValue();
      delete SFNodeVec[i];
    }

    writeFieldValues(ctx->dmSF, ctx->vStart, ctx->vEnd, sfValues.data(), "StreamFunction", ctx->viewer, step, time);
  }

  // PETSc's ExodusII VecView never calls ex_put_time (verified against PETSc 3.25.5
  // source) -- time_whole would stay at its netCDF fill value (0) for every step
  // otherwise. Patch the real time in directly through the raw ExodusII C API, using
  // the viewer's own underlying file handle.
  PetscExodusIIInt exoid;
  PetscViewerExodusIIGetId(ctx->viewer, &exoid);
  double timeValue = time;
  ex_put_time(exoid, step + 1, &timeValue);
}

void closeExodusOutput(ExodusOutputContext * ctx)
{
  PetscViewerDestroy(&ctx->viewer);
  DMDestroy(&ctx->dmU);
  DMDestroy(&ctx->dmV);
  DMDestroy(&ctx->dmP);
  DMDestroy(&ctx->dmSF);
  DMDestroy(&ctx->baseDM);
  delete ctx;
}
