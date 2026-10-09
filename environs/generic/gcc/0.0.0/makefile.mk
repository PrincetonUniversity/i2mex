# Makefile include file for generic GCC
#

ifndef MAKE_SO
	MAKE_SO = 0
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

AR = ar rcs

# Fortran flags
FFLAGS_OPT = -c -g -O2 -m64 -ffpe-trap=invalid,zero,overflow -fno-range-check -fdollar-ok $(FPIC)
FFLAGS_DBG = -c -g -O0 -m64 -ffpe-trap=invalid,zero,overflow -fno-range-check -fdollar-ok $(FPIC)
ifeq ($(shell test $(GCCVERSION) -gt 9; echo $$?),0)
        FFLAGS_OPT += -fallow-argument-mismatch -fallow-invalid-boz
        FFLAGS_DBG += -fallow-argument-mismatch -fallow-invalid-boz
endif
ifeq ($(DEBUG),0)
	FFLAGS = $(FFLAGS_OPT)
else
	FFLAGS = $(FFLAGS_DBG)
endif
MPI_FFLAGS = $(FFLAGS)

FPP = -cpp
MFLAG = -J

L_FLIBS = -L$(FLIBROOT)/lib -lgfortran -lm


# C/C++ flags
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

CPP = -cpp

L_CLIBS = -L$(FLIBROOT)/lib -lstdc++ -lgcc_s -lgcc


# For readline etc
L_EDIT = -L/usr/lib64 -lreadline -lhistory -ltermcap


# Default include directories
FINCL :=
CINCL :=


# Default definitions
FDEFS := -D__LINUX -D__UNIX -D__F90
CDEFS := -D__LINUX

FDEFS := $(FDEFS) -D__NOMDSPLUS
CDEFS := $(CDEFS) -D__NOMDSPLUS
