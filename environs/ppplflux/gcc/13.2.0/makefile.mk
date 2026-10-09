# ---------------------------------------------------------------------------
# compiler_flags.mk
#
# Standalone I2MEX compiler configuration
#
# PPPL Flux
# GCC 13.2.0
#
# Optimized for Flux AMD Zen 4 compute nodes.
#
# Environment is prepared before make with:
#
#     source environ
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
# Flux architecture optimization
#
# Flux compute nodes use AMD Zen 4 processors.
#
# -march=znver4
#     Enable instructions and code generation for AMD Zen 4.
#
# -mtune=znver4
#     Tune instruction scheduling for AMD Zen 4.
#
# Do NOT use -march=native because builds performed on a login node could
# otherwise depend on the login-node CPU rather than the Flux compute-node
# architecture.
# ===========================================================================

ARCH_FLAGS := \
	-m64 \
	-march=znver4 \
	-mtune=znver4


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
# Deliberately do not use architecture optimization beyond ARCH_FLAGS.
# -O0 makes debugging and runtime checks much easier to interpret.
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
# FLIBROOT is established by environ after loading GCC 13.2.0.
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
# Rocky Linux 9 readline/history use libtinfo.
# ===========================================================================

L_EDIT := \
	-lreadline \
	-lhistory \
	-ltinfo \
	-ldl


# ===========================================================================
# Compatibility variables used by legacy Makefiles
# ===========================================================================

FTN_MAIN :=

FTN_TRAILER := \
	$(L_FLIBS)

CPP_TRAILER := \
	$(L_CLIBS)

LOADER_TRAILER := \
	$(L_EDIT)


# ===========================================================================
# External mathematical/scientific libraries
#
# These variables are deliberately NOT redefined here.
#
# They are supplied by environ:
#
#     L_BLAS
#     L_LAPACK
#     L_NETCDF
#     L_PSPLINE
#     L_EZCDF
#     L_UFILES
#     L_MDSPLUS
#
# Their corresponding roots are also supplied by environ/module files.
# ===========================================================================
