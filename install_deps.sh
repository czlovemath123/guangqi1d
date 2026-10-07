#!/usr/bin/env bash
# Build guangqi's dependencies (OpenMPI, parallel HDF5, netlib BLAS/LAPACK,
# PETSc) into one user-owned prefix, no sudo required. Nothing outside this
# prefix is touched: the MPI is built from its source tarball too, so no
# system-wide MPI installation is needed.
#
# Usage:
#   ./install_deps.sh                 # build whatever is missing
#   ./install_deps.sh --force         # rebuild everything from scratch
#   ./install_deps.sh --prefix DIR    # install somewhere else
#   ./install_deps.sh --jobs N        # parallel make jobs (default: nproc)
#
# Tarballs: versioned names (openmpi-<ver>.tar.gz, hdf5-<ver>.tar.gz,
# petsc-<ver>.tar.gz, lapack-<ver>.tar.gz) and the version-less hdf5.tar.gz /
# petsc.tar.gz / lapack.tar.gz are both accepted from the repo root or
# $GUANGQI_SRC — the user is responsible for knowing which version their
# tarball contains. Downloading from the network is only a last resort.
#
# After it finishes: source env.sh, then make in the guangqi directory.
set -euo pipefail

OPENMPI_VERSION=${OPENMPI_VERSION:-5.0.11}
HDF5_VERSION=${HDF5_VERSION:-2.1.1}
PETSC_VERSION=${PETSC_VERSION:-3.26.0}
LAPACK_VERSION=${LAPACK_VERSION:-3.12.1}

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEPS=${GUANGQI_DEPS:-$REPO_DIR/guangqi-deps}
JOBS=$(nproc)
FORCE=0

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) DEPS=$(readlink -m "$2"); shift 2 ;;
        --force)  FORCE=1; shift ;;
        --jobs)   JOBS="$2"; shift 2 ;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "install_deps.sh: unknown option $1 (see --help)" >&2; exit 1 ;;
    esac
done

SRC_DIR=$DEPS/src
LOG_DIR=$SRC_DIR/logs
mkdir -p "$SRC_DIR" "$LOG_DIR" "$DEPS/bin" "$DEPS/lib" "$DEPS/include"

die() { echo "install_deps.sh: ERROR: $*" >&2; exit 1; }
log() { echo "==> $*"; }

# ---------------------------------------------------------------- preflight
# No MPI in the list: it is built from its tarball in the first stage.
for tool in gcc gfortran make python3 tar gzip cmake; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool not found in PATH; it is required (on Ubuntu: sudo apt install build-essential cmake python3 — the only steps that ever need root)"
done

# ------------------------------------------------------------ sources
# Look for a tarball in the repo root, then $GUANGQI_SRC; download as a last
# resort. Tarballs are copied into $SRC_DIR so builds never write to the repo.
find_tarball() {
    local name=$1 url=$2 tarball=$SRC_DIR/$1
    if [ -f "$REPO_DIR/$name" ]; then
        cp "$REPO_DIR/$name" "$tarball"
    elif [ -n "${GUANGQI_SRC:-}" ] && [ -f "$GUANGQI_SRC/$name" ]; then
        cp "$GUANGQI_SRC/$name" "$tarball"
    elif [ -f "$tarball" ]; then
        : # already fetched in a previous run
    else
        log "$name not found locally, downloading from $url"
        if command -v curl >/dev/null 2>&1; then
            curl -fL --retry 3 -o "$tarball" "$url"
        elif command -v wget >/dev/null 2>&1; then
            wget -O "$tarball" "$url"
        else
            die "$name missing and neither curl nor wget is available to download it"
        fi
    fi
}

# pick_tarball <versioned-name> <fallback-name>...: echo the first candidate
# that exists in the repo root, $GUANGQI_SRC, or $SRC_DIR; if none exists,
# echo the versioned name so find_tarball can download it. This is what makes
# the version-less hdf5.tar.gz / petsc.tar.gz / lapack.tar.gz work.
pick_tarball() {
    local name
    for name in "$@"; do
        if [ -f "$REPO_DIR/$name" ] || [ -f "$SRC_DIR/$name" ] \
           || { [ -n "${GUANGQI_SRC:-}" ] && [ -f "$GUANGQI_SRC/$name" ]; }; then
            echo "$name"; return
        fi
    done
    echo "$1"
}

# normalize_srcdir <tarball> <srcdir>: rename the extracted top-level directory
# to <srcdir> when the tarball's directory name differs from the expected
# version-named path (version-less tarballs, './'-prefixed listings).
normalize_srcdir() {
    local tarball=$1 srcdir=$2 top
    [ -d "$srcdir" ] && return 0
    top=$(tar tzf "$tarball" | grep -m1 -v '/$' \
          | awk -F/ '{for (i=1;i<=NF;i++) if ($i != "" && $i != ".") {print $i; exit}}')
    [ -n "$top" ] && [ -d "$SRC_DIR/$top" ] \
        || die "could not locate the source directory inside $tarball"
    mv "$SRC_DIR/$top" "$srcdir"
}

OPENMPI_TARBALL=openmpi-$OPENMPI_VERSION.tar.gz
HDF5_TARBALL=hdf5-$HDF5_VERSION.tar.gz
PETSC_TARBALL=petsc-$PETSC_VERSION.tar.gz
LAPACK_TARBALL=lapack-$LAPACK_VERSION.tar.gz
# Best-effort fallbacks for missing tarballs; local tarballs win.
OPENMPI_URL="https://download.open-mpi.org/release/open-mpi/v${OPENMPI_VERSION%.*}/$OPENMPI_TARBALL"
HDF5_URL="https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz"

# ---------------------------------------------------------------- helpers
# run_logged <logfile> <command...>: run a command, capture output, show the
# tail on failure.
run_logged() {
    local logfile=$1; shift
    if ! "$@" >>"$logfile" 2>&1; then
        echo "install_deps.sh: command failed: $*" >&2
        echo "--- last 40 lines of $logfile ---" >&2
        tail -n 40 "$logfile" >&2
        die "see the full log at $logfile"
    fi
}

# ---------------------------------------------------------------- OpenMPI
# Built first so that every later stage (HDF5, PETSc, and guangqi itself)
# compiles against the same MPI: hdf5.mod and petscksp.mod are only
# ABI-compatible with the mpif90 that produced them.
build_openmpi() {
    local marker=$DEPS/bin/mpif90
    if [ -x "$marker" ] && [ "$FORCE" -eq 0 ]; then
        log "OpenMPI already installed, skipping"
        return
    fi
    local tarball_name
    tarball_name=$(pick_tarball "$OPENMPI_TARBALL" openmpi.tar.gz)
    log "using tarball $tarball_name"
    find_tarball "$tarball_name" "$OPENMPI_URL"
    local srcdir=$SRC_DIR/openmpi-$OPENMPI_VERSION
    rm -rf "$srcdir"
    tar -xzf "$SRC_DIR/$tarball_name" -C "$SRC_DIR"
    normalize_srcdir "$SRC_DIR/$tarball_name" "$srcdir"
    local logfile=$LOG_DIR/openmpi-configure.log
    log "configuring OpenMPI $OPENMPI_VERSION (log: $logfile)"
    run_logged "$logfile" bash -c "cd '$srcdir' && ./configure --prefix='$DEPS'"
    logfile=$LOG_DIR/openmpi-make.log
    log "building OpenMPI with $JOBS jobs (10-20 min; log: $logfile)"
    run_logged "$logfile" bash -c "cd '$srcdir' && make -j$JOBS && make install"
    [ -x "$marker" ] || die "mpif90 not found in $DEPS/bin after install"
    log "OpenMPI $OPENMPI_VERSION installed in $DEPS"
}

# After the OpenMPI stage, $DEPS/bin must win over any system MPI so that
# HDF5 and PETSc are compiled against the same one guangqi links.
export PATH="$DEPS/bin:$PATH"

# The MPI must have Fortran bindings and wrap the same gfortran: guangqi's
# 'use mpi', hdf5.mod and petscksp.mod all come from this toolchain.
check_mpi() {
    mpif90 --version 2>&1 | grep -qiE "gfortran|GNU Fortran" \
        || die "$(command -v mpif90) does not wrap gfortran — guangqi needs an MPI built with Fortran support"
    echo "program mpi_check; use mpi; integer :: ierr; call mpi_init(ierr); call mpi_finalize(ierr); end program mpi_check" > "$SRC_DIR/mpi_check.f90"
    if ! mpif90 "$SRC_DIR/mpi_check.f90" -o "$SRC_DIR/mpi_check" 2>"$LOG_DIR/mpi_check.log"; then
        tail -n 20 "$LOG_DIR/mpi_check.log" >&2
        die "compiling a 'use mpi' test program with mpif90 failed — the MPI installation lacks usable Fortran bindings"
    fi
    # Run through mpirun (never as a singleton: some OpenMPI 5 setups hang on a
    # direct singleton MPI_Init, which guangqi never does — it always runs under mpirun).
    if ! timeout 120 mpirun --oversubscribe -np 1 "$SRC_DIR/mpi_check" >>"$LOG_DIR/mpi_check.log" 2>&1; then
        tail -n 20 "$LOG_DIR/mpi_check.log" >&2
        die "the MPI test program failed under 'mpirun -np 1' — see $LOG_DIR/mpi_check.log"
    fi
    log "MPI OK: $(mpirun --version | head -1) ($(command -v mpirun))"
}

# ---------------------------------------------------------------- HDF5
# Parallel + Fortran HDF5, built against the OpenMPI in $DEPS/bin. This is
# what produces hdf5.mod (read by guangqi's 'use hdf5') and h5dump etc.
# HDF5 >= 2.0 dropped autotools: branch to CMake there (1.x keeps configure).
build_hdf5() {
    local marker
    marker=$(ls "$DEPS"/lib/libhdf5_fortran.* 2>/dev/null | head -1 || true)
    if [ -n "$marker" ] && [ "$FORCE" -eq 0 ]; then
        log "HDF5 already installed, skipping"
        return
    fi
    # Accept a version-less tarball the user downloaded (e.g. hdf5.tar.gz);
    # the version inside is the user's responsibility.
    local tarball_name
    tarball_name=$(pick_tarball "$HDF5_TARBALL" hdf5.tar.gz)
    log "using tarball $tarball_name"
    find_tarball "$tarball_name" "$HDF5_URL"
    local srcdir=$SRC_DIR/hdf5-$HDF5_VERSION
    rm -rf "$srcdir"
    tar -xzf "$SRC_DIR/$tarball_name" -C "$SRC_DIR"
    normalize_srcdir "$SRC_DIR/$tarball_name" "$srcdir"
    local major=${HDF5_VERSION%%.*}
    if [ "$major" -ge 2 ]; then
        local logfile=$LOG_DIR/hdf5-configure.log
        log "configuring HDF5 $HDF5_VERSION with CMake + the bundled mpicc/mpif90 (log: $logfile)"
        run_logged "$logfile" cmake -S "$srcdir" -B "$srcdir/build" -G 'Unix Makefiles' \
            -DCMAKE_INSTALL_PREFIX="$DEPS" -DCMAKE_BUILD_TYPE=Release \
            -DHDF5_ENABLE_PARALLEL:BOOL=ON -DHDF5_BUILD_FORTRAN:BOOL=ON \
            -DBUILD_SHARED_LIBS:BOOL=ON -DHDF5_BUILD_HL_LIB:BOOL=ON \
            -DHDF5_BUILD_CPP_LIB:BOOL=OFF -DHDF5_BUILD_EXAMPLES:BOOL=OFF \
            -DHDF5_BUILD_UTILS:BOOL=OFF \
            -DCMAKE_C_COMPILER="$(command -v mpicc)" -DCMAKE_Fortran_COMPILER="$(command -v mpif90)"
        logfile=$LOG_DIR/hdf5-make.log
        log "building HDF5 with $JOBS jobs (log: $logfile)"
        run_logged "$logfile" cmake --build "$srcdir/build" -j "$JOBS"
        run_logged "$logfile" cmake --install "$srcdir/build"
    else
        local logfile=$LOG_DIR/hdf5-configure.log
        log "configuring HDF5 with the bundled mpicc/mpif90 (log: $logfile)"
        run_logged "$logfile" bash -c "cd '$srcdir' && CC='$(command -v mpicc)' FC='$(command -v mpif90)' ./configure --enable-parallel --enable-fortran --prefix='$DEPS'"
        logfile=$LOG_DIR/hdf5-make.log
        log "building HDF5 with $JOBS jobs (log: $logfile)"
        run_logged "$logfile" bash -c "cd '$srcdir' && make -j$JOBS && make install"
    fi
    [ -f "$DEPS/include/hdf5.mod" ] || die "hdf5.mod not found in $DEPS/include after install"
    log "HDF5 $HDF5_VERSION installed; hdf5.mod is in $DEPS/include"
}

# ---------------------------------------------------------------- LAPACK
# Netlib reference BLAS + LAPACK (the same librefblas/liblapack as in the
# user guide, built into the bundle). PETSc is configured against these, and
# guangqi links them directly too. (PETSc's own --download-f2cblaslapack does
# not compile with GCC >= 15, whose default C23 rejects its old-style
# function pointers.)
build_lapack() {
    local marker=$DEPS/lib/liblapack.a
    if [ -f "$marker" ] && [ -f "$DEPS/lib/librefblas.a" ] && [ "$FORCE" -eq 0 ]; then
        log "BLAS/LAPACK already installed, skipping"
        return
    fi
    local tarball_name
    tarball_name=$(pick_tarball "$LAPACK_TARBALL" lapack.tar.gz)
    log "using tarball $tarball_name"
    find_tarball "$tarball_name" "https://netlib.org/lapack/$LAPACK_TARBALL"
    local srcdir=$SRC_DIR/lapack-$LAPACK_VERSION
    rm -rf "$srcdir"
    tar -xzf "$SRC_DIR/$tarball_name" -C "$SRC_DIR"
    normalize_srcdir "$SRC_DIR/$tarball_name" "$srcdir"
    local logfile=$LOG_DIR/lapack-build.log
    log "building BLAS/LAPACK with gfortran (log: $logfile)"
    run_logged "$logfile" bash -c "cd '$srcdir' \
        && cp make.inc.example make.inc \
        && sed -i 's/^FFLAGS = /FFLAGS = -fPIC /; s/^FFLAGS_NOOPT = /FFLAGS_NOOPT = -fPIC /' make.inc \
        && make blaslib lapacklib -j$JOBS \
        && cp librefblas.a liblapack.a '$DEPS/lib/'"
    [ -f "$marker" ] || die "liblapack.a not found in $DEPS/lib after build"
    log "BLAS/LAPACK installed in $DEPS/lib"
}

# ---------------------------------------------------------------- PETSc
# Built against the bundled OpenMPI (via its mpicc/mpif90 wrappers) and the
# netlib BLAS/LAPACK in $DEPS/lib. --with-fortran-bindings produces the
# petscksp.mod etc. that guangqi's 'use petscksp' reads.
build_petsc() {
    local marker=$DEPS/lib/petsc/conf/variables
    if [ -f "$marker" ] && [ "$FORCE" -eq 0 ]; then
        log "PETSc already installed, skipping"
        return
    fi
    # Accept the "with docs" distribution and the version-less tarball the
    # user downloaded as well; the version inside is the user's responsibility.
    local tarball_name
    tarball_name=$(pick_tarball "$PETSC_TARBALL" petsc-with-docs-$PETSC_VERSION.tar.gz petsc.tar.gz)
    log "using tarball $tarball_name"
    find_tarball "$tarball_name" \
        "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/$PETSC_TARBALL"
    local srcdir=$SRC_DIR/petsc-$PETSC_VERSION
    rm -rf "$srcdir"
    tar -xzf "$SRC_DIR/$tarball_name" -C "$SRC_DIR"
    normalize_srcdir "$SRC_DIR/$tarball_name" "$srcdir"
    # PETSc 3.19's configure imports the stdlib module 'xdrlib', which was
    # removed in Python 3.13. Newer PETSc does not need it; only fall back to
    # an older python if one is available, otherwise run with plain python3.
    local compatpy=""
    if python3 -c "import xdrlib" 2>/dev/null; then
        compatpy=python3
    else
        local py
        for py in python3.12 python3.11 python3.10; do
            if command -v "$py" >/dev/null 2>&1 && "$py" -c "import xdrlib" 2>/dev/null; then
                compatpy=$py; break
            fi
        done
    fi
    if [ -n "$compatpy" ]; then
        mkdir -p "$SRC_DIR/shims"
        ln -sf "$(command -v "$compatpy")" "$SRC_DIR/shims/python3"
        ln -sf "$(command -v "$compatpy")" "$SRC_DIR/shims/python"
        export PATH="$SRC_DIR/shims:$PATH"
        log "using $compatpy for PETSc configure (python3 lacks xdrlib)"
    else
        log "no python with xdrlib found; using python3 as-is (fine for PETSc >= 3.20)"
    fi
    local logfile=$LOG_DIR/petsc-configure.log
    log "configuring PETSc (log: $logfile)"
    # The MPI wrappers alone tell PETSc everything about the MPI
    # installation; --with-mpi-dir must NOT be combined with them.
    # BLAS/LAPACK come from the netlib build in $DEPS/lib (see build_lapack).
    local cxxbin
    cxxbin=$(command -v mpicxx || true)
    run_logged "$logfile" bash -c "cd '$srcdir' && ./configure --prefix='$DEPS' \
        --with-cc='$(command -v mpicc)' --with-cxx='${cxxbin:-0}' --with-fc='$(command -v mpif90)' \
        --with-blas-lib='$DEPS/lib/librefblas.a' --with-lapack-lib='$DEPS/lib/liblapack.a' \
        --with-fortran-bindings=1 --with-debugging=0"
    logfile=$LOG_DIR/petsc-make.log
    log "building PETSc with $JOBS jobs (this is the long step, 30-60 min; log: $logfile)"
    run_logged "$logfile" bash -c "cd '$srcdir' && make -j$JOBS && make install"
    [ -f "$marker" ] || die "petsc conf variables not found at $marker after install"
    ls "$DEPS"/include/petscksp.mod >/dev/null 2>&1 \
        || die "petscksp.mod not found in $DEPS/include after install"
    log "PETSc $PETSC_VERSION installed (conf at $marker)"
}

build_openmpi
check_mpi
build_lapack
build_hdf5
build_petsc

log "all dependencies are in $DEPS"
log "now run:  source env.sh   (in this directory), then make"
