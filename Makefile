
# ============================================================================
# Makefile
#
# Standalone I2MEX build on Linux
#
# Supported compiler families:
#
#       GNU     : gcc / g++ / gfortran
#       Intel   : icx / icpx / ifx (or classic icc / icpc / ifort)
#       NVIDIA  : nvc / nvc++ / nvfortran
#
# Compiler-, architecture-, and system-specific configuration is supplied by:
#
#       environ
#       COMPILER_FLAGS
#
# This top-level Makefile contains no compiler-specific optimization,
# architecture, module-loading, or compute-node-specific options.
#
# Builds:
#
#       1. TRANSP-derived libraries required by I2MEX
#       2. libi2mex.a and libi2mex.so
#       3. I2MEX executable
#       4. drive executable
#       5. mex2eqs executable
#
# No TRANSP_ROOT is used.
# No installed TRANSP libraries are used.
# No MPI versions of TRANSP libraries are built here.
#
# Usage:
#
#       source environ
#       make
#
# or:
#
#       make -f Makefile
#
# Optional:
#
#       make NTASKS=16
#       make DEBUG=1
#       make show
#
# ============================================================================


# ============================================================================
# Shell
# ============================================================================

SHELL := /bin/bash


# ============================================================================
# Standalone source tree
#
# IMPORTANT:
#
# Capture the directory containing the TOP-LEVEL Makefile before including
# any other Makefiles.
#
# firstword is used deliberately: after include statements MAKEFILE_LIST
# contains additional Makefiles, but the first entry remains this top-level
# Makefile.
# ============================================================================

ROOT := $(abspath $(dir $(firstword $(MAKEFILE_LIST))))


# ============================================================================
# Required environment
# ============================================================================

ifndef LOCAL
$(error LOCAL is not defined. Run: source environ)
endif

ifndef COMPILER_FLAGS
$(error COMPILER_FLAGS is not defined. Run: source environ)
endif

ifndef UFILES_ROOT
$(error UFILES_ROOT is not defined. Run: source environ)
endif

ifndef PLASMA_STATE_ROOT
$(error PLASMA_STATE_ROOT is not defined. Run: source environ)
endif

ifndef PSPLINE_HOME
$(error PSPLINE_HOME is not defined. Run: source environ)
endif

ifndef EZCDF_HOME
$(error EZCDF_HOME is not defined. Run: source environ)
endif

ifndef MDSPLUS_DIR
$(error MDSPLUS_DIR is not defined. Run: source environ)
endif


# ============================================================================
# Compiler configuration
#
# COMPILER_FLAGS must define the compiler-dependent variables required by the
# component Makefiles, for example:
#
#       FC
#       CC
#       CXX
#       FFLAGS
#       CFLAGS
#       CXXFLAGS
#       FDEFS
#       CDEFS
#       LDFLAGS
#       L_FLIBS
#
# The selected compiler may be GNU, Intel, or NVIDIA.
# ============================================================================

include $(COMPILER_FLAGS)


# ============================================================================
# Validate compiler configuration
# ============================================================================

ifndef FC
$(error FC is not defined by environ or COMPILER_FLAGS)
endif

ifndef CC
$(error CC is not defined by environ or COMPILER_FLAGS)
endif

ifndef CXX
$(error CXX is not defined by environ or COMPILER_FLAGS)
endif


# ============================================================================
# Build options
# ============================================================================

MAKE_SO ?= 1
DEBUG   ?= 0

#
# Parallelism is allowed inside each individual library build.
#
# The outer library loop remains sequential because the order of the
# TRANSP-derived libraries is significant.
#

NTASKS ?= 8


# ============================================================================
# Source directories
# ============================================================================

I2MEX_DIR := $(ROOT)/i2mex
UTIL_DIR  := $(ROOT)/util


# ============================================================================
# TRANSP-derived libraries required by standalone I2MEX
#
# IMPORTANT:
#
# This follows the relative ordering of the TRANSP top-level Makefile.
#
# physconst is intentionally included and is built FIRST.
#
# Important later ordering:
#
#       mclib
#       rplot_mod
#       xdatmgr
#       rplot_io
#       tr_getnl
#       rp_kernel
#       trread
#       splitn
#       trxplib
#
# The outer build loop MUST remain sequential.
# ============================================================================

LIBS := \
	physconst \
	portlib \
	vaxonly \
	tokyr \
	comput \
	trdatbuf_lib \
	smlib \
	r8bloat \
	echmod_iolib \
	fluxav \
	mclib \
	rplot_mod \
	xdatmgr \
	rplot_io \
	tr_getnl \
	rp_kernel \
	trread \
	splitn \
	trxplib


# ============================================================================
# I2MEX executables
#
# i2mex and drive both use util/drive.f90 as their main program.
# ============================================================================

EXES := \
	i2mex \
	drive \
	mex2eqs


# ============================================================================
# Local I2MEX link dependencies
#
# This is LINK order, not BUILD order.
#
# All TRANSP-derived libraries are taken from:
#
#       $(LOCAL)/lib
#
# physconst is intentionally included.
# ============================================================================

L_I2MEX_LOCAL_DEP := \
	-L$(LOCAL)/lib \
	-ltrxplib \
	-ltrread \
	-lr8bloat \
	-lrp_kernel \
	-ltr_getnl \
	-lsplitn \
	-lechmod_iolib \
	-ltrdatbuf_lib \
	-lrplot_io \
	-lrplot_mod \
	-lxdatmgr \
	-lfluxav \
	-lmclib \
	-ltokyr \
	-lvaxonly \
	-lsmlib \
	-lcomput \
	-lportlib \
	-lphysconst


# ============================================================================
# External dependencies
#
# These are supplied by environ and/or COMPILER_FLAGS and are not built here.
#
# They may point to compiler-specific installations appropriate for the
# selected GNU, Intel, or NVIDIA toolchain.
# ============================================================================

L_I2MEX_EXTERNAL_DEP := \
	$(L_PLASMA_STATE) \
	$(L_UFILES) \
	$(L_MDSPLUS) \
	$(L_PSPLINE) \
	$(L_EZCDF) \
	$(L_NETCDF) \
	$(L_LAPACK) \
	$(L_BLAS) \
	$(L_FLIBS)


L_I2MEX_DEP := \
	$(L_I2MEX_LOCAL_DEP) \
	$(L_I2MEX_EXTERNAL_DEP)


# ============================================================================
# Linux linker search paths
#
# -Wl passes the following option to the system linker.
#
# These options are accepted by the normal Linux compiler drivers used by:
#
#       GNU
#       Intel
#       NVIDIA HPC SDK
#
# rpath:
#       stored in the resulting ELF executable and used at runtime.
#
# rpath-link:
#       used by the linker to locate indirect shared-library dependencies
#       while linking; it is not a runtime search path.
#
# The compiler-specific runtime-library paths, if required, belong in
# L_FLIBS / LDFLAGS from COMPILER_FLAGS rather than here.
# ============================================================================

L_I2MEX_RPATH := \
	-Wl,-rpath,$(LOCAL)/lib \
	-Wl,-rpath-link,$(LOCAL)/lib


# ============================================================================
# Utility object directory
# ============================================================================

UTIL_OBJDIR := $(LOCAL)/obj/i2mex_util

DRIVE_OBJ   := $(UTIL_OBJDIR)/drive.o
MEX2EQS_OBJ := $(UTIL_OBJDIR)/mex2eqs.o


# ============================================================================
# Utility include/module paths
#
# Fortran .mod files are compiler-specific. Therefore all external packages
# referenced here must have been built with a compiler compatible with FC.
# ============================================================================

UTIL_FINCL := \
	-I$(LOCAL)/mod \
	-I$(PLASMA_STATE_ROOT)/mod \
	-I$(PSPLINE_HOME)/include \
	-I$(EZCDF_HOME)/mod


# ============================================================================
# Export variables required by component Makefiles
# ============================================================================

export LOCAL
export COMPILER_FLAGS

export MAKE_SO
export DEBUG

export L_I2MEX_LOCAL_DEP
export L_I2MEX_EXTERNAL_DEP
export L_I2MEX_DEP
export L_I2MEX_RPATH


# ============================================================================
# Targets
# ============================================================================

.PHONY: \
	all \
	check \
	libs \
	i2mex \
	exes \
	show \
	clean \
	realclean \
	clobber \
	depend


# ============================================================================
# Normal build sequence
#
#       check
#         |
#        libs
#         |
#     libi2mex
#         |
#    executables
#
# There is intentionally NO separate TRANSP MPI-library build stage.
# ============================================================================

all: exes


# ============================================================================
# Show configuration
# ============================================================================

show:
	@echo
	@echo "======================================================================"
	@echo " Standalone I2MEX build configuration"
	@echo "======================================================================"
	@echo
	@echo "ROOT              = $(ROOT)"
	@echo "I2MEX_DIR         = $(I2MEX_DIR)"
	@echo "UTIL_DIR          = $(UTIL_DIR)"
	@echo
	@echo "LOCAL             = $(LOCAL)"
	@echo "COMPILER_FLAGS    = $(COMPILER_FLAGS)"
	@echo
	@echo "MAKE_SO           = $(MAKE_SO)"
	@echo "DEBUG             = $(DEBUG)"
	@echo "NTASKS            = $(NTASKS)"
	@echo
	@echo "FC                = $(FC)"
	@echo "CC                = $(CC)"
	@echo "CXX               = $(CXX)"
	@echo
	@echo "FFLAGS            = $(FFLAGS)"
	@echo "CFLAGS            = $(CFLAGS)"
	@echo "CXXFLAGS          = $(CXXFLAGS)"
	@echo "LDFLAGS           = $(LDFLAGS)"
	@echo
	@echo "UFILES_ROOT       = $(UFILES_ROOT)"
	@echo "PLASMA_STATE_ROOT = $(PLASMA_STATE_ROOT)"
	@echo "PSPLINE_HOME      = $(PSPLINE_HOME)"
	@echo "EZCDF_HOME        = $(EZCDF_HOME)"
	@echo "MDSPLUS_DIR       = $(MDSPLUS_DIR)"
	@echo
	@echo "Libraries:"
	@for l in $(LIBS); do \
		echo "  $$l"; \
	done
	@echo
	@echo "Executables:"
	@for e in $(EXES); do \
		echo "  $$e"; \
	done
	@echo
	@echo "L_I2MEX_LOCAL_DEP:"
	@echo "  $(L_I2MEX_LOCAL_DEP)"
	@echo
	@echo "L_I2MEX_EXTERNAL_DEP:"
	@echo "  $(L_I2MEX_EXTERNAL_DEP)"
	@echo
	@echo "L_I2MEX_RPATH:"
	@echo "  $(L_I2MEX_RPATH)"
	@echo
	@echo "======================================================================"


# ============================================================================
# Check/create build directories
# ============================================================================

check:
	@mkdir -p $(LOCAL)/bin
	@mkdir -p $(LOCAL)/exe
	@mkdir -p $(LOCAL)/lib
	@mkdir -p $(LOCAL)/mod
	@mkdir -p $(LOCAL)/obj
	@mkdir -p $(LOCAL)/tmp
	@mkdir -p $(LOCAL)/work
	@mkdir -p $(UTIL_OBJDIR)

	@set -e; \
	for l in $(LIBS); do \
		if [ ! -d "$(ROOT)/$${l}" ]; then \
			echo "ERROR: missing source directory: $(ROOT)/$${l}"; \
			exit 1; \
		fi; \
		if [ ! -f "$(ROOT)/$${l}/$${l}_makefile.mk" ]; then \
			echo "ERROR: missing Makefile: $(ROOT)/$${l}/$${l}_makefile.mk"; \
			exit 1; \
		fi; \
	done

	@if [ ! -d "$(I2MEX_DIR)" ]; then \
		echo "ERROR: missing I2MEX source directory: $(I2MEX_DIR)"; \
		exit 1; \
	fi

	@if [ ! -f "$(I2MEX_DIR)/i2mex_makefile.mk" ]; then \
		echo "ERROR: missing Makefile: $(I2MEX_DIR)/i2mex_makefile.mk"; \
		exit 1; \
	fi

	@if [ ! -d "$(UTIL_DIR)" ]; then \
		echo "ERROR: missing utility source directory: $(UTIL_DIR)"; \
		exit 1; \
	fi

	@if [ ! -f "$(UTIL_DIR)/drive.f90" ]; then \
		echo "ERROR: missing source: $(UTIL_DIR)/drive.f90"; \
		exit 1; \
	fi

	@if [ ! -f "$(UTIL_DIR)/mex2eqs.f90" ]; then \
		echo "ERROR: missing source: $(UTIL_DIR)/mex2eqs.f90"; \
		exit 1; \
	fi

	@set -e; \
	for l in $(LIBS); do \
		cd "$(ROOT)/$${l}"; \
		$(MAKE) --no-print-directory \
			-f $${l}_makefile.mk \
			MAKE_SO=$(MAKE_SO) \
			DEBUG=$(DEBUG) \
			check; \
	done

	@cd "$(I2MEX_DIR)" && \
		$(MAKE) --no-print-directory \
			-f i2mex_makefile.mk \
			MAKE_SO=$(MAKE_SO) \
			DEBUG=$(DEBUG) \
			check


# ============================================================================
# Build TRANSP-derived libraries
#
# IMPORTANT:
#
# The OUTER loop is sequential.
#
# Individual library Makefiles may build their own source files in parallel,
# but the next library does not begin until the current library is complete.
# ============================================================================

libs: check
	@set -e; \
	for l in $(LIBS); do \
		echo; \
		echo "============================================================"; \
		echo " Building $${l}"; \
		echo "============================================================"; \
		cd "$(ROOT)/$${l}"; \
		$(MAKE) --no-print-directory \
			-f $${l}_makefile.mk \
			-j $(NTASKS) \
			MAKE_SO=$(MAKE_SO) \
			DEBUG=$(DEBUG) \
			$${l}; \
	done


# ============================================================================
# Build libi2mex
# ============================================================================

i2mex: libs
	@echo
	@echo "============================================================"
	@echo " Building I2MEX library"
	@echo "============================================================"
	@cd "$(I2MEX_DIR)" && \
		$(MAKE) --no-print-directory \
			-f i2mex_makefile.mk \
			-j $(NTASKS) \
			MAKE_SO=$(MAKE_SO) \
			DEBUG=$(DEBUG) \
			i2mex


# ============================================================================
# Compile drive.f90
#
# This object is used for BOTH:
#
#       $(LOCAL)/bin/i2mex
#       $(LOCAL)/bin/drive
# ============================================================================

$(DRIVE_OBJ): $(UTIL_DIR)/drive.f90
	@mkdir -p $(UTIL_OBJDIR)
	@echo
	@echo "============================================================"
	@echo " Compiling drive.f90"
	@echo "============================================================"
	$(FC) \
		$(FFLAGS) \
		$(FDEFS) \
		$(UTIL_FINCL) \
		$< \
		-o $@


# ============================================================================
# Compile mex2eqs.f90
# ============================================================================

$(MEX2EQS_OBJ): $(UTIL_DIR)/mex2eqs.f90
	@mkdir -p $(UTIL_OBJDIR)
	@echo
	@echo "============================================================"
	@echo " Compiling mex2eqs.f90"
	@echo "============================================================"
	$(FC) \
		$(FFLAGS) \
		$(FDEFS) \
		$(UTIL_FINCL) \
		$< \
		-o $@


# ============================================================================
# Link I2MEX executable
#
# drive.f90 supplies the main program.
# ============================================================================

$(LOCAL)/bin/i2mex: $(DRIVE_OBJ)
	@mkdir -p $(LOCAL)/bin
	@echo
	@echo "============================================================"
	@echo " Linking executable: i2mex"
	@echo "============================================================"
	$(FC) \
		-o $@ \
		$(DRIVE_OBJ) \
		-L$(LOCAL)/lib \
		-li2mex \
		$(L_I2MEX_DEP) \
		$(L_I2MEX_RPATH)


# ============================================================================
# Link drive
# ============================================================================

$(LOCAL)/bin/drive: $(DRIVE_OBJ)
	@mkdir -p $(LOCAL)/bin
	@echo
	@echo "============================================================"
	@echo " Linking executable: drive"
	@echo "============================================================"
	$(FC) \
		-o $@ \
		$(DRIVE_OBJ) \
		-L$(LOCAL)/lib \
		-li2mex \
		$(L_I2MEX_DEP) \
		$(L_I2MEX_RPATH)


# ============================================================================
# Link mex2eqs
# ============================================================================

$(LOCAL)/bin/mex2eqs: $(MEX2EQS_OBJ)
	@mkdir -p $(LOCAL)/bin
	@echo
	@echo "============================================================"
	@echo " Linking executable: mex2eqs"
	@echo "============================================================"
	$(FC) \
		-o $@ \
		$(MEX2EQS_OBJ) \
		-L$(LOCAL)/lib \
		-li2mex \
		$(L_I2MEX_DEP) \
		$(L_SGLIB) \
		$(L_I2MEX_RPATH)


# ============================================================================
# Build executables
# ============================================================================

exes: i2mex \
	$(LOCAL)/bin/i2mex \
	$(LOCAL)/bin/drive \
	$(LOCAL)/bin/mex2eqs
	@echo
	@echo "============================================================"
	@echo " Standalone I2MEX build complete"
	@echo "============================================================"
	@echo
	@echo "Libraries:"
	@echo "  $(LOCAL)/lib/libi2mex.a"
	@echo "  $(LOCAL)/lib/libi2mex.so"
	@echo
	@echo "Executables:"
	@echo "  $(LOCAL)/bin/i2mex"
	@echo "  $(LOCAL)/bin/drive"
	@echo "  $(LOCAL)/bin/mex2eqs"
	@echo


# ============================================================================
# Clean
# ============================================================================

clean:
	@set -e; \
	for l in $(LIBS); do \
		echo "Cleaning $${l}"; \
		mkdir -p "$(LOCAL)/mod"; \
		if [ -d "$(ROOT)/$${l}" ] && \
		   [ -f "$(ROOT)/$${l}/$${l}_makefile.mk" ]; then \
			cd "$(ROOT)/$${l}"; \
			$(MAKE) --no-print-directory \
				-f $${l}_makefile.mk \
				MAKE_SO=$(MAKE_SO) \
				DEBUG=$(DEBUG) \
				clean; \
		fi; \
	done

	@if [ -d "$(I2MEX_DIR)" ] && \
	    [ -f "$(I2MEX_DIR)/i2mex_makefile.mk" ]; then \
		echo "Cleaning i2mex"; \
		mkdir -p "$(LOCAL)/mod"; \
		cd "$(I2MEX_DIR)" && \
		$(MAKE) --no-print-directory \
			-f i2mex_makefile.mk \
			MAKE_SO=$(MAKE_SO) \
			DEBUG=$(DEBUG) \
			clean; \
	fi

	@echo "Cleaning I2MEX utility objects/executables"
	@rm -rf $(UTIL_OBJDIR)
	@rm -f $(LOCAL)/bin/i2mex
	@rm -f $(LOCAL)/bin/drive
	@rm -f $(LOCAL)/bin/mex2eqs


# ============================================================================
# Real clean
# ============================================================================

realclean:
	@set -e; \
	for l in $(LIBS); do \
		echo "Realclean $${l}"; \
		mkdir -p "$(LOCAL)/mod"; \
		if [ -d "$(ROOT)/$${l}" ] && \
		   [ -f "$(ROOT)/$${l}/$${l}_makefile.mk" ]; then \
			cd "$(ROOT)/$${l}"; \
			$(MAKE) --no-print-directory \
				-f $${l}_makefile.mk \
				MAKE_SO=$(MAKE_SO) \
				DEBUG=$(DEBUG) \
				realclean; \
		fi; \
	done

	@if [ -d "$(I2MEX_DIR)" ] && \
	    [ -f "$(I2MEX_DIR)/i2mex_makefile.mk" ]; then \
		echo "Realclean i2mex"; \
		mkdir -p "$(LOCAL)/mod"; \
		cd "$(I2MEX_DIR)" && \
		$(MAKE) --no-print-directory \
			-f i2mex_makefile.mk \
			MAKE_SO=$(MAKE_SO) \
			DEBUG=$(DEBUG) \
			realclean; \
	fi

	@rm -rf $(LOCAL)/bin
	@rm -rf $(LOCAL)/exe
	@rm -rf $(LOCAL)/lib
	@rm -rf $(LOCAL)/mod
	@rm -rf $(LOCAL)/obj
	@rm -rf $(LOCAL)/tmp
	@rm -rf $(LOCAL)/work


# ============================================================================
# Clobber standalone build products
# ============================================================================

clobber:
	@echo "Removing standalone I2MEX build products under:"
	@echo "  $(LOCAL)"
	@rm -rf $(LOCAL)/bin
	@rm -rf $(LOCAL)/exe
	@rm -rf $(LOCAL)/lib
	@rm -rf $(LOCAL)/mod
	@rm -rf $(LOCAL)/obj
	@rm -rf $(LOCAL)/tmp
	@rm -rf $(LOCAL)/work


# ============================================================================
# Dependency generation
# ============================================================================

depend:
	@set -e; \
	for l in $(LIBS); do \
		echo "Checking dependency target for $${l}"; \
		cd "$(ROOT)/$${l}"; \
		if $(MAKE) --no-print-directory \
			-f $${l}_makefile.mk \
			-n depend >/dev/null 2>&1; then \
			echo "Generating dependencies for $${l}"; \
			$(MAKE) --no-print-directory \
				-f $${l}_makefile.mk \
				MAKE_SO=$(MAKE_SO) \
				DEBUG=$(DEBUG) \
				depend; \
		else \
			echo "No depend target for $${l}; skipping"; \
		fi; \
	done
