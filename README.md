# Guangqi 1D — AT 2025abao light-curve reproduction

Guangqi (光启) is a radiation-hydrodynamics code named after Xu Guangqi
(徐光启), the 16th-century Chinese mathematician and astronomer who
collaborated with Matteo Ricci. It is originally a **2D code** published as
*[Guangqi: A two-dimensional radiation hydrodynamic code with realistic
equation of states](https://arxiv.org/abs/2602.12545)* (Chen & Bai,
arXiv:2602.12545); this repository hosts the **1D (spherically symmetric)
version**, which specializes in radiation transport with realistic
H/He equations of state.

The `abao.py` pipeline in `modules/lrne/` reproduces the light-curve fits of
the AT 2025abao luminous red nova study,
[arXiv:2610.05769](https://arxiv.org/abs/2610.05769): it runs two complete
simulations of the ejecta (a low-speed and a high-speed realization), and
compares their bolometric light curves against the observed AT 2025abao UVOIR
data.

## 1. Dependencies

The Fortran code needs:

| dependency | version used here | notes |
|------------|-------------------|-------|
| gcc / gfortran | any recent (tested with GCC 13) | Fortran 2003+, preprocessor |
| OpenMPI | **5.0.11** | with Fortran bindings; **built from source by `install_deps.sh`** — no system MPI needed |
| HDF5 | **2.1.1** | parallel, with Fortran bindings |
| PETSc | **3.26.0** | with Fortran bindings (implicit radiation solver) |
| BLAS / LAPACK | **3.12.1** (netlib reference) | |
| make, cmake, python3, tar, gzip | system versions | cmake is required by HDF5 2.x; no root access needed |

The Python driver (`modules/lrne/abao.py`) needs **Python 3 with numpy and
matplotlib** (any virtualenv works; e.g. `/home/zhuo/git/myenv/bin/python` on
the author's machine).

All four libraries (OpenMPI, HDF5, PETSc, LAPACK/BLAS) are built from source
into the repository — you do **not** need to install any of them system-wide.

## 2. What the shell scripts do

- **`install_deps.sh`** builds a complete local dependency folder
  `guangqi-deps/` inside the repository. It looks for the source tarballs
  (`openmpi-5.0.11.tar.gz`, `hdf5-2.1.1.tar.gz`, `petsc-3.26.0.tar.gz`,
  `lapack-3.12.1.tar.gz`; version-less names like `hdf5.tar.gz` also work) in
  the repository root or in `$GUANGQI_SRC`, downloads them only as a last
  resort, and builds OpenMPI first, then netlib BLAS/LAPACK, parallel HDF5,
  and PETSc against that exact MPI — so everything is ABI-consistent. It is
  idempotent (already-built stages are skipped; `--force` rebuilds), logs each
  step to `guangqi-deps/src/logs/`, and never touches anything outside the
  prefix. Options: `--prefix DIR`, `--jobs N`, `--force`.
- **`env.sh`** puts the bundle on your environment: it prepends
  `guangqi-deps/bin` (mpiexec, mpif90, h5fc, …) to `PATH`,
  `guangqi-deps/lib` to `LD_LIBRARY_PATH`, and exports `PETSC_DIR`. Source it
  in every shell before building or running.

## 3. Running `abao.py` from scratch

```bash
# 0) get the source tarballs into this repository (or set GUANGQI_SRC
#    to the directory holding them)
cp openmpi-5.0.11.tar.gz hdf5-2.1.1.tar.gz petsc-3.26.0.tar.gz \
   lapack-3.12.1.tar.gz  /path/to/guangqi1d

cd /path/to/guangqi1d

# 1) build the local dependency folder (OpenMPI + LAPACK/BLAS + HDF5 + PETSc)
./install_deps.sh --jobs 24          # ~30-60 min, once per checkout

# 2) activate the environment
source env.sh

# 3) run the campaign (builds the code automatically on first use)
cd modules/lrne
python abao.py all                   # init -> run -> bc -> lc -> rhol -> vhist -> escmass
```

`abao.py run` points the `modules/problem` symlink at `lrne`, runs
`make clean && make` if the `guangqi` binary is missing (use
`python abao.py run --rebuild` to force a rebuild), and executes each case
with `mpiexec -np 1` — **the lrne module runs on a single MPI rank only**.
Each 160-day simulation takes a few minutes on one core. The steps can also
be run individually (`init`, `run`, `bc`, `lc`, `rhol`, `vhist`, `escmass`).

Results, per case (`modules/lrne/lowspeed/` and `modules/lrne/highspeed/`):

- `out/lrneNNNNN.h5` (+ `.xdmf`) — snapshot frames for ParaView; the `out/`
  directory is created automatically by the code if missing.
- `history.data` — time series (column 0: time [s], column 3: outer
  luminosity [erg/s]).
- `run.log`, `last.dat` — simulation log and last frame index.

Figures, in `modules/lrne/`: `best_lc_comparison.png` (model light curves
overplotted on the observed AT 2025abao UVOIR data — the main result of
[arXiv:2610.05769](https://arxiv.org/abs/2610.05769)), `best_boundcond.png`
(the two ejecta tables), `pictures/rhoL.png` (Kippenhahn-style evolution
maps), `best_vhist.png` and `escape_mass.png`.

## 4. Further documentation

- [USER_GUIDE.md](USER_GUIDE.md) — code-wide guide: prerequisites, build
  system, compile-time switches, input decks, running.
- [modules/lrne/README.md](modules/lrne/README.md) — the LRNe physics module:
  model description, input decks, outputs, and the full `abao.py` reference.
