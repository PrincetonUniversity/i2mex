# Makefile include file for PPPL Intel 2019
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
FFLAGS_OPT = -c -g -O2 -nowarn -ftz -traceback -align dcommons -heap-arrays 20 $(FPIC)
FFLAGS_DBG = -c -g -O0 -w -check bounds -check pointer -check uninit -fpe0 -traceback -align dcommons -heap-arrays 20 $(FPIC)
ifeq ($(DEBUG),0)
	FFLAGS = $(FFLAGS_OPT)
else
	FFLAGS = $(FFLAGS_DBG)
endif

FPP = -fpp
MFLAG = -module

L_FLIBS = -L$(FLIBROOT)/lib -lifport -limf -lifcore -lm -lirc


# C/C++ flags
CFLAGS_OPT = -c -g -O2 -w -mp $(FPIC)
CFLAGS_DBG = -c -g -O0 -w -mp $(FPIC)
CXXFLAGS_OPT = -c -g -O2 -w -mp $(FPIC)
CXXFLAGS_DBG = -c -g -O0 -w -mp $(FPIC)
ifeq ($(DEBUG),0)
	CFLAGS = $(CFLAGS_OPT)
	CXXFLAGS = $(CXXFLAGS_OPT)
else
	CFLAGS = $(CFLAGS_DBG)
	CXXFLAGS = $(CXXFLAGS_DBG)
endif

CPP = -cpp

L_CLIBS = -L/usr/lib64 -lc -lstdc++


# For readline etc
L_EDIT = -L/usr/lib64 -lreadline -lhistory -ltermcap


# Default include directories
FINCL :=
CINCL :=


# Default definitions
FDEFS := -D__LINUX -D__UNIX -D__F90
CDEFS := -D__LINUX

# MDS+ not installed on stellar
FDEFS := $(FDEFS) -D__NOMDSPLUS
CDEFS := $(CDEFS) -D__NOMDSPLUS
