#!/usr/bin/env bash
# Source this file before building or running guangqi:
#   source env.sh
# It makes the OpenMPI/HDF5/PETSc/BLAS-LAPACK built by ./install_deps.sh
# findable. No system-wide MPI is required: mpirun/mpif90 come from the
# bundle itself.
_GUANGQI_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export GUANGQI_DEPS=${GUANGQI_DEPS:-$_GUANGQI_REPO/guangqi-deps}

if [ ! -x "$GUANGQI_DEPS/bin/mpirun" ]; then
    echo "env.sh: no mpirun in $GUANGQI_DEPS/bin — run ./install_deps.sh first" >&2
    return 1
fi
if [ ! -e "$GUANGQI_DEPS/lib/libhdf5.so" ] && [ ! -e "$GUANGQI_DEPS/lib/libhdf5.a" ]; then
    echo "env.sh: no HDF5 in $GUANGQI_DEPS — run ./install_deps.sh first" >&2
    return 1
fi

case ":$PATH:" in
    *":$GUANGQI_DEPS/bin:"*) ;;
    *) export PATH="$GUANGQI_DEPS/bin:$PATH" ;;
esac
case ":${LD_LIBRARY_PATH:-}:" in
    *":$GUANGQI_DEPS/lib:"*) ;;
    *) export LD_LIBRARY_PATH="$GUANGQI_DEPS/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
esac
export PETSC_DIR=$GUANGQI_DEPS

echo "env.sh: using guangqi-deps at $GUANGQI_DEPS (MPI: $(command -v mpirun))"
