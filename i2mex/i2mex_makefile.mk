# Makefile to build libi2mex.a and libi2mex.so for the TRANSP build system
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := i2mex

OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR  := $(LOCAL)/lib/lib$(NAME).a
LIBSO  := $(LOCAL)/lib/lib$(NAME).so


# ----------------------------------------------------------------------
# Sources and modules
# ----------------------------------------------------------------------

MODS := cont_mod.mod freeqbe_mod.mod i2mex_mod.mod imex_ode_mod.mod

SRCM := $(MODS:%.mod=%.f90)
SRCF := $(SRCM) $(filter-out $(SRCM),$(wildcard *.f90))
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)

SRCC := $(wildcard *.c)
OBJC := $(SRCC:%.c=$(OBJDIR)/%.o)


# ----------------------------------------------------------------------
# Compiler flags
# ----------------------------------------------------------------------

FDEF := $(FDEFS)

FINC := $(FINCL) \
	-I$(NETCDF_FORTRAN_HOME)/include \
	-I$(PSPLINE_HOME)/include \
	-I$(EZCDF_HOME)/mod \
	-I$(PLASMA_STATE_ROOT)/mod

FFLG := $(FFLAGS)

ifeq ($(MAKE_SO),1)
	FFLG += -fPIC
endif


CDEF := $(CDEFS)
CINC := $(CINCL)
CFLG := $(CFLAGS)

ifeq ($(MAKE_SO),1)
	CFLG += -fPIC
endif


# ----------------------------------------------------------------------
# Shared-library dependencies
#
# These are TRANSP libraries used directly by I2MEX.
# They must be recorded as dependencies of libi2mex.so when MAKE_SO=1.
# ----------------------------------------------------------------------

L_I2MEX_DEP := \
	-L$(LOCAL)/lib \
	-ltrxplib \
	-ltrread \
	-lr8bloat \
	-lsmlib \
	-lfluxav \
	-lmclib \
	-lcomput \
	-lvaxonly \
	-lportlib


# ----------------------------------------------------------------------
# External dependencies
# ----------------------------------------------------------------------

L_I2MEX_DEP += \
	$(L_PLASMA_STATE) \
	$(L_MDSPLUS) \
	$(L_PSPLINE) \
	$(L_EZCDF) \
	$(L_NETCDF) \
	$(L_BLAS) \
	$(L_FLIBS)


# ----------------------------------------------------------------------
# Runtime search paths
#
# $ORIGIN allows libi2mex.so to find TRANSP shared libraries installed
# in the same $(LOCAL)/lib directory.
# ----------------------------------------------------------------------

RPATH_I2MEX := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath-link,$(LOCAL)/lib


# ----------------------------------------------------------------------
# Targets
# ----------------------------------------------------------------------

.PHONY: all clean realclean clobber check depend

all: check $(NAME)


check:
	@mkdir -p $(LOCAL)/lib
	@mkdir -p $(LOCAL)/obj/$(NAME)
	@mkdir -p $(LOCAL)/mod


# ----------------------------------------------------------------------
# Build static and optional shared libraries
# ----------------------------------------------------------------------

$(NAME): $(OBJF) $(OBJC)
	@echo "Building $(NAME) static library"
	@rm -f $(LIBAR)
	@$(AR) $(LIBAR) $(OBJF) $(OBJC)
ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@rm -f $(LIBSO)
	@$(FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) $(OBJC) \
		$(L_I2MEX_DEP) \
		$(RPATH_I2MEX)
endif


# ----------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------

include $(NAME)_depend.mk


# ----------------------------------------------------------------------
# Compilation
# ----------------------------------------------------------------------

$(OBJDIR)/%.o: %.f90 | check
	@$(FC) $(FFLG) \
		$< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mod


$(OBJDIR)/%.o: %.c | check
	@$(CC) $(CFLG) \
		$< \
		-o $@ \
		$(CDEF) \
		$(CINC)


# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod && rm -f $(MODS)


realclean: clean
	@rm -f $(LIBAR)
	@rm -f $(LIBSO)


clobber: realclean


# ----------------------------------------------------------------------
# Generate Fortran dependencies
# ----------------------------------------------------------------------

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
