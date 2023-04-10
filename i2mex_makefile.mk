# Makefile to build libi2mex.a for the TRANSP build system
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := i2mex
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := cont_mod.mod freeqbe_mod.mod i2mex_mod.mod
SRCM := $(MODS:%.mod=%.f90)
SRCF := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)

SRCC := $(wildcard *.c)
OBJC := $(SRCC:%.c=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL) -I$(NETCDF_FORTRAN_HOME)/include -I$(PSPLINE_HOME)/include -I$(EZCDF_HOME)/mod
FFLG := $(FFLAGS)

CDEF := $(CDEFS)
CINC := $(CINCL)
CFLG := $(CFLAGS)


.PHONY: clean realclean clobber check

all:	check $(NAME)

check:
	@test -d $(LOCAL)/lib || mkdir -p $(LOCAL)/lib
	@test -d $(LOCAL)/obj/$(NAME) || mkdir -p $(LOCAL)/obj/$(NAME)

$(NAME): $(OBJF) $(OBJC)
	@echo "Building $(NAME) static library"
	@$(AR) $(LIBAR) $(OBJF) $(OBJC)
ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@$(FC) -shared $(LDFLAGS) -o $(LIBSO) $(OBJF) $(OBJC)
endif

include $(NAME)_depend.mk
$(OBJDIR)/%.o: %.f90
	@$(FC)  $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

$(OBJDIR)/%.o: %.c
	@$(CC)  $(CFLG) $< -o $@ $(CDEF) $(CINC)

clean: 
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod ; rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
