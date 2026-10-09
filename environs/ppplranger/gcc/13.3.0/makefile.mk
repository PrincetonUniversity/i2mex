# ---------------------------------------------------------------------------
# compiler_flags.mk
#
# Standalone I2MEX compiler configuration
#
# PPPL Ranger
# GCC 13.3.0
#
# Optimized for Ranger Intel Xeon compute nodes.
#
# Environment is prepared before make with:
#
#     source environs/ppplranger/gcc/13.3.0/environ
#
# No TRANSP installation is required.
# ---------------------------------------------------------------------------


# ===========================================================================
# Build options
# ===========================================================================

ifndef MAKE_SO
	MAKE_SO := 1
endif

ifndef DEBUG
	DEBUG := 0
endif


# ===========================================================================
# Compilers
#
# Normally exported by environ:
#
#     FC=gfortran
#     CC=gcc
#     CXX=g++
# ===========================================================================

FC  ?= gfortran
CC  ?= gcc
CXX ?= g++


# ===========================================================================
# Archiver
# ===========================================================================

AR := ar rcs


# ===========================================================================
# Shared-library support
# ===========================================================================

FPIC :=

ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
endif

ifeq ($(MAKE_SO),TRUE)
	FPIC := -fPIC
endif

ifeq ($(MAKE_SO),true)
	FPIC := -fPIC
endif


# ===========================================================================
# Ranger architecture optimization
#
# Ranger compute nodes use Intel Xeon processors from the Skylake/Cascade
# Lake server generation.
#
# GCC uses the skylake-avx512 target for this processor family.
#
# -march=skylake-avx512
#     Generate code for Intel Skylake server-class processors with AVX-512.
#     This target is also appropriate for Cascade Lake processors.
#
# -mtune=skylake-avx512
#     Tune instruction scheduling and code generation for this architecture.
#
# Do NOT use:
#
#     -march=native
#
# because a build performed on a login node could otherwise be optimized for
# the login-node processor instead of the Ranger compute-node architecture.
#
# If Ranger later contains compute nodes with an older/different CPU
# architecture, ARCH_FLAGS can be changed here without modifying the
# standalone I2MEX Makefile.
# ===========================================================================

ARCH_FLAGS := \
	-m64 \
	-march=skylake-avx512 \
	-mtune=skylake-avx512


# ===========================================================================
# Fortran preprocessing / module output
# ===========================================================================

FPP   := -cpp
MFLAG := -J


# ===========================================================================
# Fortran flags
#
# -c
#     Compile only.  The individual library Makefiles perform linking.
#
# -fallow-argument-mismatch
#     Required for legacy Fortran interfaces where argument types/ranks
#     are not fully consistent.
#
# -fallow-invalid-boz
#     Permit legacy BOZ literal usage accepted by older Fortran compilers.
#
# -fno-range-check
#     Required by legacy source containing constants that modern gfortran
#     would otherwise reject at compile time.
#
# -fdollar-ok
#     Permit '$' in legacy Fortran identifiers.
#
# -fPIC
#     Added when MAKE_SO=1 for shared-library construction.
# ===========================================================================

FFLAGS_COMMON := \
	-c \
	$(ARCH_FLAGS) \
	-fallow-argument-mismatch \
	-fallow-invalid-boz \
	-fno-range-check \
	-fdollar-ok \
	$(FPIC)


# ---------------------------------------------------------------------------
# Optimized production build
# ---------------------------------------------------------------------------

FFLAGS_OPT := \
	$(FFLAGS_COMMON) \
	-O2 \
	-g


# ---------------------------------------------------------------------------
# Debug build
#
# Keep the same instruction-set target as the production build, but disable
# compiler optimization so debugging and runtime diagnostics remain useful.
# ---------------------------------------------------------------------------

FFLAGS_DBG := \
	$(FFLAGS_COMMON) \
	-O0 \
	-g \
	-ffpe-trap=invalid,zero,overflow \
	-fcheck=all \
	-fbacktrace


ifeq ($(DEBUG),0)
	FFLAGS := $(FFLAGS_OPT)
else
	FFLAGS := $(FFLAGS_DBG)
endif


# ===========================================================================
# C flags
# ===========================================================================

CFLAGS_COMMON := \
	-c \
	$(ARCH_FLAGS) \
	$(FPIC)


CFLAGS_OPT := \
	$(CFLAGS_COMMON) \
	-O2 \
	-g


CFLAGS_DBG := \
	$(CFLAGS_COMMON) \
	-O0 \
	-g


ifeq ($(DEBUG),0)
	CFLAGS := $(CFLAGS_OPT)
else
	CFLAGS := $(CFLAGS_DBG)
endif


# ===========================================================================
# C++ flags
# ===========================================================================

CXXFLAGS_COMMON := \
	-c \
	$(ARCH_FLAGS) \
	$(FPIC)


CXXFLAGS_OPT := \
	$(CXXFLAGS_COMMON) \
	-O2 \
	-g


CXXFLAGS_DBG := \
	$(CXXFLAGS_COMMON) \
	-O0 \
	-g


ifeq ($(DEBUG),0)
	CXXFLAGS := $(CXXFLAGS_OPT)
else
	CXXFLAGS := $(CXXFLAGS_DBG)
endif


# ===========================================================================
# Preprocessor aliases used by legacy Makefiles
# ===========================================================================

CPP := -cpp


# ===========================================================================
# Preprocessor definitions
# ===========================================================================

FDEFS := \
	-D__LINUX \
	-D__UNIX \
	-D__F90

CDEFS := \
	-D__LINUX


# ===========================================================================
# Default include flags
#
# Individual library Makefiles add their own include paths.
#
# NetCDF Fortran modules on Ranger are installed under:
#
#     /usr/lib64/gfortran/modules
#
# The individual Makefiles may use NETCDF_MODDIR directly when required.
# ===========================================================================

FINCL :=
CINCL :=


# ===========================================================================
# General linker flags
# ===========================================================================

LDFLAGS :=


# ===========================================================================
# Fortran runtime
#
# FLIBROOT is established by environ after enabling GCC Toolset 13.
#
# For GCC 13.3.0:
#
#     FLIBROOT=$(dirname $(gfortran -print-file-name=libgfortran.so))
#
# This ensures that the same GCC runtime used for compilation is found when
# the standalone I2MEX shared libraries and executables are loaded.
# ===========================================================================

L_FLIBS := \
	-Wl,-rpath,$(FLIBROOT) \
	-L$(FLIBROOT) \
	-lgfortran \
	-lm


# ===========================================================================
# C/C++ runtime libraries
# ===========================================================================

L_CLIBS := \
	-Wl,-rpath,$(FLIBROOT) \
	-L$(FLIBROOT) \
	-lstdc++ \
	-lgcc_s \
	-lgcc


# ===========================================================================
# Readline / terminal libraries
#
# Ranger provides these libraries under /usr/lib64.
#
# Ranger uses libtermcap for the terminal compatibility library in the
# existing TRANSP/Plasma-State environment.
# ===========================================================================

L_EDIT := \
	-L/usr/lib64 \
	-lreadline \
	-lhistory \
	-ltermcap \
	-ldl


# ===========================================================================
# Compatibility variables used by legacy TRANSP-derived Makefiles
# ===========================================================================

FTN_MAIN :=

FTN_TRAILER := \
	$(L_FLIBS)

CPP_TRAILER := \
	$(L_CLIBS)

LOADER_TRAILER := \
	$(L_EDIT)


# ===========================================================================
# NetCDF / HDF5 defaults
#
# Ranger RPM installations use /usr for NetCDF and HDF5 and place the
# gfortran module files in /usr/lib64/gfortran/modules.
#
# These are defaults only.  Values exported by environ take precedence.
# ===========================================================================

NETCDF_MODDIR ?= /usr/lib64/gfortran/modules
NETCDF_F_HOME ?= /usr
NETCDF_C_HOME ?= /usr
HDF5_HOME     ?= /usr


# ===========================================================================
# External mathematical/scientific libraries
#
# These variables are deliberately NOT redefined here.
#
# They are supplied by the Ranger environ:
#
#     L_BLAS
#     L_LAPACK
#     L_NETCDF
#     L_PSPLINE
#     L_EZCDF
#     L_UFILES
#     L_MDSPLUS
#     L_PLASMA_STATE
#
# Their corresponding roots are also supplied by environ:
#
#     PSPLINE_HOME
#     EZCDF_HOME
#     UFILES_ROOT
#     PLASMA_STATE_ROOT
#     MDSPLUS_HOME
#     NETCDF_F_HOME
#     NETCDF_C_HOME
#     HDF5_HOME
#
# The compiler configuration therefore remains independent of the installed
# locations of these external packages.
# ===========================================================================
