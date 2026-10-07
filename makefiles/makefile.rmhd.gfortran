All: guangqi
OBJ_DIR = obj
VPATH = src:modules/problem

objects = $(addprefix $(OBJ_DIR)/, datastructure.o tree_module.o mathlib.o phylib.o	\
	constant_gamma_eos.o eos_h2_HI.o eos_h_he.o eos_h2_HI_HII.o eos_HI_HII.o eos_analytic.o eos.o)
objects += $(addprefix $(OBJ_DIR)/, problem.o)
objects += $(addprefix $(OBJ_DIR)/, communication_1d.o communication.o load_module.o	\
	passivescalars.o radiation_common_functions.o boundary.o amr_module.o io_out.o io_in.o	\
	hllc.o eos_hllc_analytic.o hydro.o gravity.o petsc_fld_1d.o petsc_fld.o \
	radiation.o source_control.o recon_evolve.o muscl.o hydroscheme.o rmhd.o)
ieos ?= 1
ieosmodule ?=1
iopacity ?= 1
isolver ?= 1
irecord ?= 0
imodify ?= 0
usersource ?= 0
isource_order ?= 0
user_amr ?= 0
include modules/problem/makefile.problem
ompfflag = -fopenmp
guangqi_flags = $(foreach var,\
	ieos ieosmodule iopacity isolver openmp irecord imodify isource_order usersource user_amr,\
	-D$(var)=$($(var)))

# Dependencies built by ./install_deps.sh into <repo>/guangqi-deps (override:
#   make DEPS_DIR=/path/to/deps, or export GUANGQI_DEPS before sourcing env.sh).
GUANGQI_DEPS ?= $(CURDIR)/guangqi-deps
DEPS_DIR = $(GUANGQI_DEPS)
PETSC_DIR = $(DEPS_DIR)
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
FC = mpif90

#FFLAGS += ${PETSC_FC_INCLUDES} -cpp -ffree-line-length-512 -fcheck=all -g -O2 $(guangqi_flags)
#FFLAGS += ${PETSC_FC_INCLUDES} -cpp -ffree-line-length-512 -fcheck=all -g -fbacktrace -O0 $(guangqi_flags)
FFLAGS += ${PETSC_FC_INCLUDES} -cpp -ffree-line-length-512 -O3	$(guangqi_flags)

#guangqi's environment
HDF5 = $(DEPS_DIR)
openmpi = $(DEPS_DIR)
include_path += -I$(HDF5)/include -I$(openmpi)/include -I$(openmpi)/lib -I$(PETSC_DIR)/include
blas = -lrefblas
lapack = -llapack
LIB_DIR = -L$(HDF5)/lib
LIBS = -lhdf5 -lhdf5_fortran $(lapack) $(blas) -lm

COMPILE = $(FC) -fPIC $(FFLAGS) $(include_path) -J$(OBJ_DIR) -I$(OBJ_DIR)
LINK = $(FC) -fPIC $(FFLAGS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)
	$(COMPILE) -c $< -o $@
guangqi: $(objects)
	${LINK} $(objects) $(LIB_DIR) $(LIBS) ${PETSC_LIB} -o guangqi

.PHONY: clean
clean::
	rm -f $(objects) $(OBJ_DIR)/*.mod *.mod guangqi gmon.out analysis.txt
