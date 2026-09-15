# Guangqi 1D User Guide

This guide describes how to compile and run the 1D version of Guangqi, a
radiation-hydrodynamics code that solves the radiation transport equations
coupled to hydrodynamics with complex equations of state.

## Repository layout

```
├── Makefile                 # top-level build file (builds the `guangqi` executable)
├── rmhd.f90                 # main program
├── src/                     # solver source files (hydro, radiation, eos, io, ...)
├── makefiles/               # machine-specific makefile examples
├── modules/                 # problem modules (initial/boundary conditions)
│   ├── lrne/                #   luminous red nova (ejecta of CEE) setup
│   ├── acc_planet/          #   accreting gas giant setup
│   └── problem -> lrne      #   symlink selecting the ACTIVE module
├── tables/                  # opacity tables (required at runtime)
└── scripts/                 # python analysis / table utilities
```

A *problem module* defines the physical problem: `problem.f90` (initial and
boundary conditions), `makefile.problem` (compile-time switches), input decks
(`global.data`, `problem.data`), and any module data files. The code builds one
module at a time through the `modules/problem` symlink.

## Prerequisites

- **gfortran** with **OpenMPI** (`mpif90`, `mpiexec`)
- **PETSc** built with Fortran support (used for the implicit radiation solver)
- **HDF5** with Fortran bindings (parallel output)
- **LAPACK / BLAS**
- **Python 3** for the helper/analysis scripts (numpy, matplotlib)

The `Makefile` contains machine paths that you likely need to set. They
default to placeholders and can be overridden via environment variables or on
the make command line:

```make
PETSC_DIR ?= /path/to/petsc        # your PETSc installation (required)
HDF5      ?= /usr/local/hdf5       # your HDF5 installation
openmpi   ?= /usr/local/openmpi    # your OpenMPI installation
```

so e.g. `export PETSC_DIR=$HOME/opt/petsc` once, and both `make` and
`abao.py run` (which honors `PETSC_DIR`, `HDF5_DIR`, `OPENMPI_DIR`) pick it
up. `make` fails with a clear error if `PETSC_DIR` does not point at a PETSc
installation.

`makefiles/makefile.rmhd.gfortran` shows a complete example of a gfortran
configuration.

## Compile-time switches

Each module ships a `makefile.problem` that pins the compile-time physics
switches, e.g. `modules/lrne/makefile.problem`:

```make
ieos = 2            # equation of state selection
ieosmodule = 3      # EOS module variant
isolver = 2         # hydrodynamic solver
iopacity = 1        # opacity treatment
irecord = 1         # write history.data time series
isource_order = 0   # source-term integration order
```

Available variables: `ieos ieosmodule iopacity isolver irecord imodify
isource_order usersource user_amr`. They are passed to the compiler as
`-D` flags; you can override any of them on the command line, e.g.
`make isolver=1`. Switching modules changes these flags automatically via the
module's `makefile.problem`.

## Building

```bash
mkdir -p obj        # the object directory is not tracked by git
make                # builds the executable ./guangqi
make clean          # removes objects, .mod files and the executable
```

To work on a different problem, first repoint the module symlink and rebuild
from clean:

```bash
ln -sfn acc_planet modules/problem
make clean && make
```

## Running

The code reads its input decks from the **current working directory** and
writes all output there, so you run the executable from inside the module
directory.

1. Prepare the run directory (shown for `lrne`):

   ```bash
   cd modules/lrne
   cp output_var_info.dat.empty output_var_info.dat   # list of variables stored in the frames
   ```

   `output_var_info.dat` is required at startup: its first line is the number
   of recorded variables, followed by one variable name per line. Edit it to
   control what goes into the output frames.

2. Edit `global.data` (the global input deck, Fortran namelist format):

   - `$meshinfo` — domain (`n_domain`), unit lengths (`lengthscale`,
     `timescale`), final time (`tfinal`), `CFL`, boundary condition types,
     number of saved frames (`nframe`), refinement settings, and restart
     options.
   - `$phyinfo` — composition and physics switches.
   - `$global_parameters` — gravity, geometry (`igeometry`: 0 Cartesian,
     1 cylindrical, 2 polar), cooling, radiation scheme (`iradiation`:
     0 none, 1 grey dust, 2 short characteristics, 4 FLD), resolution
     (`nx`, `ny`, `blk_size_nx/ny`), and `path_tables`.

   > **Important:** `path_tables` is an absolute path into the repository's
   > `tables/opacity` directory. Update it to match your checkout, e.g.
   > `path_tables='/home/you/guangqi1d/tables/opacity'`.

   Boundary condition types (used by `hydro_bound_type`, `rad_bound_type`,
   `passive_bound_type`): `1` transmissive, `2` reflective, `3` periodic,
   `4` diode, `8`/`9` specified (time-dependent inner boundaries; see
   `problem.f90` for details).

3. Edit `problem.data` — the module-specific namelist (`&parameters_1d`)
   with the physical setup of the problem (masses, floors, boundary
   configuration, recording cadence, ...). See the comments in the file and
   in the module's `problem.f90`.

4. Launch with MPI:

   ```bash
   mpiexec -np 4 ../../guangqi        # or copy ./guangqi here and run ./guangqi
   ```

   The number of MPI ranks is flexible; blocks are distributed statically at
   start-up.

### Output

- `out/NNNNN.h5` — HDF5 snapshot frames written every `tfinal/nframe`, with
  `.xdmf` companions for ParaView/Visit.
- `history.data` — time series of diagnostics (when `irecord=1`).
- `last.dat` — index of the latest written frame.

### Restarting

In `global.data` set `restart = .true.` and `restart_iframe` to the frame
number to resume from (the corresponding `.h5` must exist in `out/`), then
rerun with the same executable and number of frames adjusted through
`nframe`.

## The acc_planet module

`modules/acc_planet` studies an accreting gas giant and includes a small
batch driver, `run.sh`, which runs a scan over model parameters:

```bash
cd modules/acc_planet
ln -sfn acc_planet modules/problem     # from the repo root, then rebuild
./run.sh input.txt
```

Each non-empty line of `input.txt` contains four numbers:

```
acc_rate m_planet r_planet kh_timescale
```

For every line the driver rewrites the corresponding entries of
`problem.data` (`acc_rate` in 10⁻³ M⊕/yr, `m_planet` in M_J, `kh_timescale`
in 10⁵ yr) via `modify_parameters.py`, runs `mpiexec -np 4 ./guangqi`, and
archives the outputs with `rename.py`.

## Post-processing

Python utilities live in `scripts/` (`eos.py`, `phy_const.py`,
`table_gen.py`, `assemble_1d_data.py`, ...) and in each module (e.g.
`modules/lrne/kippenhahn2.py` for Kippenhahn diagrams). The HDF5 frames are
best visualized with ParaView through the generated `.xdmf` files.
