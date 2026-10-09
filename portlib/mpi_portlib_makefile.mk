# Makefile to build libmpi_portlib.a and libmpi_portlib.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := mpi_portlib
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so


# ======================================================================
# Sources and objects
# ======================================================================

MODS := mpi_proc_data.mod mpi_env_mod.mod logmod.mod execsystem.mod mmpi.mod

SRCM := $(MODS:%.mod=%.f90)

SF90 := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)

SRCF := $(SF90) $(SFOR)

OBJF := $(SF90:%.f90=$(OBJDIR)/%.o)
OBJF += $(SFOR:%.for=$(OBJDIR)/%.o)

SRCC := $(wildcard *.c)
OBJC := $(SRCC:%.c=$(OBJDIR)/%.o)


# ======================================================================
# Compiler flags
# ======================================================================

FDEF := $(FDEFS) -D__MPI

FINC := $(FINCL) \
	-I$(LOCAL)/mpi_mod \
	-I$(LOCAL)/mod

CDEF := $(CDEFS) \
	-D__MPI \
	-D__GETLINE_EDITOR

CINC := $(CINCL) \
	-I../incl_cpp

ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(MPI_FFLAGS) $(FPIC)
CFLG := $(MPI_CFLAGS) $(FPIC)


# ======================================================================
# Shared-library dependencies
# ======================================================================

L_MPI_PORTLIB_DEP := \
	-L$(PLASMA_STATE_ROOT)/lib \
	-lps_portlib \
	-lreadline


# ======================================================================
# Runtime search paths
# ======================================================================

RPATH_MPI_PORTLIB := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(PLASMA_STATE_ROOT)/lib


# ======================================================================
# Targets
# ======================================================================

.PHONY: all check clean realclean clobber depend

all: check $(NAME)


# ======================================================================
# Directory setup
# ======================================================================

check:
	@mkdir -p $(LOCAL)/lib
	@mkdir -p $(OBJDIR)
	@mkdir -p $(LOCAL)/mpi_mod


# ======================================================================
# Build static and shared libraries
# ======================================================================

$(NAME): $(OBJF) $(OBJC)
	@echo "Building $(NAME) static library"
	@$(AR) $(LIBAR) $(OBJF) $(OBJC)

ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@$(MPI_FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) $(OBJC) \
		$(L_MPI_PORTLIB_DEP) \
		$(RPATH_MPI_PORTLIB)
endif


# ======================================================================
# Fortran dependencies
# ======================================================================

include $(NAME)_depend.mk


# ======================================================================
# Compilation rules
# ======================================================================

$(OBJDIR)/%.o: %.for
	@$(MPI_FC) $(FPP) $(FFLG) $< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mpi_mod

$(OBJDIR)/%.o: %.f90
	@$(MPI_FC) $(FPP) $(FFLG) $< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mpi_mod

$(OBJDIR)/%.o: %.c
	@$(MPI_CC) -cpp $(CFLG) $< \
		-o $@ \
		$(CDEF) \
		$(CINC)


# ======================================================================
# Clean
# ======================================================================

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mpi_mod && rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean


# ======================================================================
# Generate dependencies
# ======================================================================

depend:
	@makedepf90 -b OBJDIR $(FDEF) $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
