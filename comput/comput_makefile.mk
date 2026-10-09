# Makefile to build libcomput.a and libcomput.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME   := comput
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR  := $(LOCAL)/lib/lib$(NAME).a
LIBSO  := $(LOCAL)/lib/lib$(NAME).so


# ---------------------------------------------------------------------------
# Fortran modules
# ---------------------------------------------------------------------------

MNMS := \
	catfus_mod.mod \
	lincir_contour_mod.mod \
	namerefs.mod \
	periodic_table_mod.mod \
	plasma_hash_code.mod \
	r8bsmoo_mod.mod \
	submsg_mod.mod \
	ode_mod.mod

MODS := \
	catfus_mod.o \
	lincir_contour.o \
	namerefs.o \
	periodic_table.o \
	plasma_hash_code.o \
	r8bsmoo_mod.o \
	submsg.o \
	ode.o


# ---------------------------------------------------------------------------
# Sources and objects
# ---------------------------------------------------------------------------

SRCM := $(MODS:%.o=%.f90)

SF90 := $(SRCM) $(filter-out $(SRCM),$(wildcard *.f90))
SFOR := $(wildcard *.for)

SRCF := $(SF90) $(SFOR)

OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

SRCC := $(wildcard *.c)
OBJC := $(SRCC:%.c=$(OBJDIR)/%.o)


# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------

FDEF := $(FDEFS)
FINC := $(FINCL)
FFLG := $(FFLAGS)

CDEF := $(CDEFS)
CINC := $(CINCL) -I../incl_cpp
CFLG := $(CFLAGS)


# ---------------------------------------------------------------------------
# Shared-library dependencies
#
# libcomput.so directly requires libportlib.so.
#
# $ORIGIN makes the runtime loader search the directory containing
# libcomput.so itself.  In the standalone build both libcomput.so and
# libportlib.so are installed in:
#
#     $(LOCAL)/lib
#
# Therefore libcomput.so does not depend on LD_LIBRARY_PATH to locate
# libportlib.so.
#
# libps_portlib.so is NOT linked here explicitly.  It is brought in through
# the dependency chain of libportlib.so.
# ---------------------------------------------------------------------------

L_COMPUT_DEP := \
	-L$(LOCAL)/lib \
	-lportlib

RPATH_COMPUT := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath-link,$(LOCAL)/lib


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

.PHONY: all clean realclean clobber check depend

all: check $(NAME)


# ---------------------------------------------------------------------------
# Create output directories
# ---------------------------------------------------------------------------

check:
	@test -d $(LOCAL)/lib || mkdir -p $(LOCAL)/lib
	@test -d $(LOCAL)/mod || mkdir -p $(LOCAL)/mod
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)


# ---------------------------------------------------------------------------
# Build libraries
# ---------------------------------------------------------------------------

$(NAME): $(OBJF) $(OBJC)
	@echo "Building $(NAME) static library"
	@$(AR) $(LIBAR) $(OBJF) $(OBJC)

ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@$(FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) $(OBJC) \
		$(L_COMPUT_DEP) \
		$(RPATH_COMPUT)
endif


# ---------------------------------------------------------------------------
# Dependency file
# ---------------------------------------------------------------------------

include $(NAME)_depend.mk


# ---------------------------------------------------------------------------
# Fortran compilation
# ---------------------------------------------------------------------------

$(OBJDIR)/%.o: %.for
	@$(FC) \
		$(FFLG) \
		$< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mod

$(OBJDIR)/%.o: %.f90
	@$(FC) \
		$(FPP) \
		$(FFLG) \
		$< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mod


# ---------------------------------------------------------------------------
# C compilation
# ---------------------------------------------------------------------------

$(OBJDIR)/%.o: %.c
	@$(CC) \
		$(CFLG) \
		$< \
		-o $@ \
		$(CDEF) \
		$(CINC)


# ---------------------------------------------------------------------------
# Cleaning
# ---------------------------------------------------------------------------

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod && rm -f $(MNMS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean


# ---------------------------------------------------------------------------
# Generate Fortran dependency file
# ---------------------------------------------------------------------------

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
