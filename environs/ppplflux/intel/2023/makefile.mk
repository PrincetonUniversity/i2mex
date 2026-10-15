
# ======================================================================
# I2MEX compiler flags
#
# PPPL Flux - Rocky Linux 9
# Intel oneAPI 2023.2.0 Classic compilers
#
# Compilers:
#   ifort
#   icc
#   icpc
#
# Supports:
#   Static libraries
#   Shared libraries (MAKE_SO=1)
#   Debug builds (DEBUG=1)
# ======================================================================

# ======================================================================
# Build configuration
# ======================================================================

MAKE_SO ?= 0
DEBUG ?= 0

FPIC :=

ifneq ($(filter 1 TRUE true,$(MAKE_SO)),)
    FPIC := -fPIC
endif

# ======================================================================
# Compiler selection
# ======================================================================

FC = ifort
FC90 = $(FC)

CC = icc
CXX = icpc

MPI_FC = $(if $(strip $(MPIF90)),$(MPIF90),mpif90)
MPI_CC = $(if $(strip $(MPICC)),$(MPICC),mpicc)
MPI_CXX = $(if $(strip $(MPICXX)),$(MPICXX),mpicxx)

# ======================================================================
# Archive utility
#
# Preserve compatibility with recipes using:
#   $(AR) $@ $(OBJF)
# ======================================================================

AR = ar rcs

# ======================================================================
# Floating-point options
#
# -fp-model precise : preserve more predictable FP behavior
# -no-ftz           : preserve subnormal floating-point values
#
# Do not use -xHost for binaries shared across Flux nodes.
# ======================================================================

FLUX_FLG = -fp-model precise -no-ftz

# ======================================================================
# Memory model
#
# Retained from original I2MEX configuration.
#
# Remove -mcmodel=medium if not needed for large static data.
# ======================================================================

MEMMODEL = -mcmodel=medium

# ======================================================================
# Fortran compiler flags
# ======================================================================

FFLAGS_OPT = \
    -c \
    -O2 \
    -nowarn \
    -align dcommons \
    -heap-arrays 20 \
    $(MEMMODEL) \
    $(FLUX_FLG) \
    $(FPIC)

FFLAGS_DBG = \
    -c \
    -g \
    -O0 \
    -check bounds \
    -check pointer \
    -check uninit \
    -fpe0 \
    -traceback \
    -align dcommons \
    -heap-arrays 20 \
    -init=zero \
    $(MEMMODEL) \
    $(FLUX_FLG) \
    $(FPIC)

ifeq ($(DEBUG),0)
    FFLAGS = $(FFLAGS_OPT)
else
    FFLAGS = $(FFLAGS_DBG)
endif

MPI_FFLAGS = $(FFLAGS)

# ======================================================================
# Fortran preprocessing and module output
# ======================================================================

FPP = -fpp
MFLAG = -module

# ======================================================================
# Intel Fortran runtime libraries
#
# FLIBROOT already ends in:
#   compiler/lib/intel64_lin
# ======================================================================

L_FLIBS = \
    -L$(FLIBROOT) \
    -lifport \
    -lintlc \
    -limf \
    -lifcore \
    -lirc \
    -lsvml \
    -lm \
    -ldl

# ======================================================================
# C compiler flags
# ======================================================================

CFLAGS_OPT = \
    -c \
    -O2 \
    -fp-model precise \
    $(FPIC)

CFLAGS_DBG = \
    -c \
    -g \
    -O0 \
    -fp-model precise \
    $(FPIC)

# ======================================================================
# C++ compiler flags
# ======================================================================

CXXFLAGS_OPT = \
    -c \
    -O2 \
    -fp-model precise \
    $(FPIC)

CXXFLAGS_DBG = \
    -c \
    -g \
    -O0 \
    -fp-model precise \
    $(FPIC)

ifeq ($(DEBUG),0)
    CFLAGS = $(CFLAGS_OPT)
    CXXFLAGS = $(CXXFLAGS_OPT)
else
    CFLAGS = $(CFLAGS_DBG)
    CXXFLAGS = $(CXXFLAGS_DBG)
endif

MPI_CFLAGS = $(CFLAGS)
MPI_CXXFLAGS = $(CXXFLAGS)

# ======================================================================
# C preprocessing
#
# Preserve existing variable for compatibility.
# ======================================================================

CPP = -cpp

# ======================================================================
# C/C++ runtime libraries
# ======================================================================

L_CLIBS = -lstdc++ -lm -ldl

# ======================================================================
# Readline and terminal libraries
# ======================================================================

L_EDIT = -lreadline -lhistory -ltinfo

# ======================================================================
# Include directories
# ======================================================================

FINCL :=
CINCL :=

# ======================================================================
# Preprocessor definitions
# ======================================================================

FDEFS := -D__LINUX -D__UNIX -D__F90
CDEFS := -D__LINUX

# ======================================================================
# Optional definitions
# ======================================================================

# FDEFS += -D__NOMDSPLUS
# CDEFS += -D__NOMDSPLUS
