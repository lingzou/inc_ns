# inc_ns

Cell-centered finite-volume solver for 2D incompressible Navier-Stokes flow on
unstructured triangular meshes, following:

> Amaresh Dalal, V. Eswaran, and G. Biswas, "A Finite-Volume Method for Navier-Stokes
> Equations on Unstructured Meshes," *Numerical Heat Transfer, Part B*, 54: 238-259, 2008.

## Dependencies

### Fluent2DMesh

Mesh I/O (reading Fluent/Gambit `.msh` files) lives in a separate sibling repo,
[Fluent2DMesh](../Fluent2DMesh), expected to be checked out next to this one:

```
projects/
├── inc_ns/          (this repo)
└── Fluent2DMesh/
```

`CMakeLists.txt` finds it via the `FLUENT2DMESH_DIR` cache variable, which defaults to
`../Fluent2DMesh` relative to this repo -- no editing needed for the side-by-side layout
above. If you check it out somewhere else, point at it explicitly:

```sh
cmake -S . -B build -DFLUENT2DMESH_DIR=<path-to-Fluent2DMesh> ...
```

### PETSc, built with ExodusII support

The solver uses PETSc for its linear solves (KSP) and for output (the `exodusii`
`PetscViewer`). You need to build/install PETSc yourself; this repo does not vendor or
build one. Two things matter about that install:

- It must have ExodusII support, i.e. configured with `--download-exodusii` (or a system
  ExodusII it was pointed at) -- this is *not* the default for many prebuilt/conda PETSc
  packages, since `--download-exodusii` is opt-in. Without it, `PetscViewerExodusIIOpen`
  fails at runtime with something like `"Unknown PetscViewer type given: exodusii"`.
- `PETSC_DIR`/`PETSC_ARCH` need to point at it when you configure this project (either
  as environment variables or `-DPETSC_DIR=...`/`-DPETSC_ARCH=...` on the `cmake` line).

## Build

```sh
cmake -S . -B build -DPETSC_DIR=<path-to-petsc> -DPETSC_ARCH=<petsc-arch>
cmake --build build
```

`compile.sh` is a working example pinned to this machine's own PETSc install -- treat it
as a template, not something to run verbatim elsewhere.

## Run

```sh
./build/inc_exec -i input/cylinder.i
```

Everything is driven by the `-i <file>.i` input file (see **Input file structure**
below). Output is written to `output/`, named from the input file's own basename (e.g.
`cylinder.i` -> `cylinder_out.*`):

- `<problem>_out.e` -- Exodus output (open in ParaView).
- `<problem>_out_<step>.vtu` -- one VTU per time step, kept alongside Exodus for now.
- `<problem>_out_convergence.csv` -- per-time-step `u`/`v` velocity L2-error history.

Regression tests (`exodiff` against `expected/*_out.e`, via SEACAS's `exodiff`):

```sh
./runTests.sh
```

## GUI

```sh
python gui/gui.py
```

A small Tkinter app: `File -> Open` loads an input file into a read-only preview, and
`Run` executes it and live-plots the `u`/`v` L2 error convergence history (tailing the
CSV above) as the simulation progresses. Needs `tkinter` and `matplotlib`; no PETSc or
build step required to just run the GUI itself.

## Input file structure

A case is entirely defined by its `-i <file>.i` input file -- mesh, timestepping, and
every zone's boundary condition -- parsed with PETSc's own options-file mechanism (plain
`-key value` lines, `#` comments allowed). See `input/cylinder.i` and `input/cavity.i`
for complete, working examples. Required top-level keys:

```
-problem    <name>            # a label, for your own reference
-mesh       <path/to/mesh>    # resolved relative to the directory inc_exec runs from
-dt         <time step>
-visc       <kinematic viscosity>
-n_timestep <number of time steps>
```

Every face zone in the mesh (`inc_exec` prints each as `Face Zone <id>` on startup) must
also have a `-bc_<zone>_type <type>` entry, plus whatever values that type needs:

| type                   | values                                        | meaning |
|------------------------|------------------------------------------------|---------|
| `interior`             | none                                           | shared by two cells; not a boundary |
| `velocity_inlet`       | `-bc_<zone>_u`, `-bc_<zone>_v` (required)      | both velocity components fixed |
| `wall`                 | `-bc_<zone>_u`, `-bc_<zone>_v` (optional, default 0) | no-slip; a *moving* wall (e.g. a lid-driven cavity's lid) is just a `wall` with a non-zero value |
| `horizontal_slip_wall` | none                                           | `v = 0` only, `u` left free |
| `vertical_slip_wall`   | none                                           | `u = 0` only, `v` left free |
| `pressure_outlet`      | `-bc_<zone>_p` (required)                      | fixed pressure |

If no zone provides a pressure B.C. at all (e.g. the lid-driven cavity, walls all the way
around), the pressure equation is only defined up to a constant; the solver detects this
automatically and anchors one cell's pressure to zero, no input-file setting needed.

A missing `-bc_<zone>_type`, an unrecognized type name, or a missing required value for a
given type makes `inc_exec` exit immediately with an error naming the specific
zone/option at fault, rather than running with silently wrong physics.
