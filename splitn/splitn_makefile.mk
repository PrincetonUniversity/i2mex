# Makefile to build libsplitn.a
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := splitn
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := splitn_module.mod
SRCM := $(MODS:%.mod=%.f90)
SF90 := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)
SRCF := $(SF90) $(SFOR)
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)
ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(FFLAGS) $(FPIC)

L_SPLITN_DEP := \
	-L$(LOCAL)/lib \
	-lportlib \
	-L$(PLASMA_STATE_ROOT)/lib \
	-lps_portlib

RPATH_SPLITN := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(PLASMA_STATE_ROOT)/lib

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
	@$(FC) -shared $(LDFLAGS) -o $(LIBSO) $(OBJF) \
	$(L_SPLITN_DEP) \
	$(RPATH_SPLITN)
endif

include $(NAME)_depend.mk
$(OBJDIR)/%.o: %.for
	@$(FC) $(FPP) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod
$(OBJDIR)/%.o: %.f90
	@$(FC) $(FPP) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod ; rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
