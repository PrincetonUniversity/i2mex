# Makefile to build librplot_mod.a
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := rplot_mod
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := cplotr_mod.mod extrac2_mod.mod mfblok_mod.mod nltrdat_mod.mod plfmpa_mod.mod rpcalc_mod.mod
SRCM := $(MODS:%.mod=%.f90)
SRCF := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)
ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(FFLAGS) $(FPIC)



.PHONY: clean realclean clobber check

all:	check $(NAME)

check:
	@test -d $(LOCAL)/lib || mkdir -p $(LOCAL)/lib
	@test -d $(LOCAL)/obj/$(NAME) || mkdir -p $(LOCAL)/obj/$(NAME)

$(NAME): $(OBJF)
	@echo "Building $(NAME) static library"
	@$(AR) $(LIBAR) $(OBJF)
ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@$(FC) -shared $(LDFLAGS) -o $(LIBSO) $(OBJF)
endif

include $(NAME)_depend.mk
$(OBJDIR)/%.o: %.f90
	@$(FC) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod -I$(NETCDF_FORTRAN_HOME)/include

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod ; rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
