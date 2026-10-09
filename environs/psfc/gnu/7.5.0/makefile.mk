# Makefile include file for PSFC GCC 7.5.0
#

ifndef MAKE_SO
	MAKE_SO = 0
endif
ifndef OMP_ACTIVATE
	OMP_ACTIVATE = 0
endif
ifndef DEBUG
	DEBUG = 0
endif

ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
endif
ifeq ($(MAKE_SO),TRUE)
	FPIC := -fPIC
endif

AR = ar -rU

# Fortran flags
FC = mpif90
FC90 = $(FC)
MPI_FC = mpif90

FFLAGS_OPT = -c -g -O2 -m64 -ffpe-trap=invalid,zero,overflow -fno-range-check -fdollar-ok $(FPIC)
FFLAGS_DBG = -c -g -O0 -m64 -ffpe-trap=invalid,zero,overflow -fno-range-check -fdollar-ok $(FPIC)
ifeq ($(DEBUG),0)
	FFLAGS = $(FFLAGS_OPT)
else
	FFLAGS = $(FFLAGS_DBG)
endif
MPI_FFLAGS = $(FFLAGS)

FOPENMP = -fopenmp

FPP = -cpp
MFLAG = -J

L_FLIBS = -L$(FLIBROOT)/lib -lgfortran -lm


# C/C++ flags
CC = mpicc
CXX= mpicxx
MPI_CC = mpicc
MPI_CXX = mpicxx

CFLAGS_OPT = -c -g -O2 -m64 $(FPIC)
CFLAGS_DBG = -c -g -O0 -m64 $(FPIC)
CXXFLAGS_OPT = -c -g -O2 -m64 $(FPIC)
CXXFLAGS_DBG = -c -g -O0 -m64 $(FPIC)
ifeq ($(DEBUG),0)
	CFLAGS = $(CFLAGS_OPT)
	CXXFLAGS = $(CXXFLAGS_OPT)
else
	CFLAGS = $(CFLAGS_DBG)
	CXXFLAGS = $(CXXFLAGS_DBG)
endif
MPI_CFLAGS = $(CFLAGS)
MPI_CXXFLAGS = $(CXXFLAGS)

L_CLIBS = -L$(FLIBROOT)/lib -lstdc++ -lgcc_s -lgcc


# For readline etc
L_EDIT = -L/usr/lib64 -lreadline -lhistory -ltermcap


# Default include directories
FINCL :=
CINCL :=


# Default definitions
FDEFS := -D__LINUX -D__UNIX -D__F90
CDEFS := -D__LINUX

FDEFS := $(FDEFS) -D__NOCQL3D -D__NOGENRAY -D__NOTGLF -D__NOGLF23 -D__NONEO -D__NOMMM
CDEFS := $(FDEFS) -D__NOCQL3D -D__NOGENRAY -D__NOTGLF -D__NOGLF23 -D__NONEO -D__NOMMM

NO_TGLF = 1
NO_GLF23 = 1
NO_NEO = 1
NO_MMM = 1
NO_CQL3D = 1
NO_GENRAY = 1
