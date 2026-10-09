# Makefile to build libechmod_iolib.a
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := echmod_iolib
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

SF90 := $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)
SRCF := $(SF90) $(SFOR)
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL) -I$(EZCDF_HOME)/mod -I$(PLASMA_STATE_ROOT)/mod
ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(FFLAGS) $(FPIC)

L_ECHMOD_IOLIB_DEP := \
	-L$(LOCAL)/lib \
	-lportlib \
	-lcomput \
	$(L_PLASMA_STATE) \
	$(L_EZCDF) \
	-L$(UFILES_ROOT)/lib \
	-lmdstransp

RPATH_ECHMOD_IOLIB := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(PLASMA_STATE_ROOT)/lib \
	-Wl,-rpath,$(EZCDF_HOME)/lib \
	-Wl,-rpath,$(UFILES_ROOT)/lib

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
	$(L_ECHMOD_IOLIB_DEP) \
	$(RPATH_ECHMOD_IOLIB)
endif

include $(NAME)_depend.mk
$(OBJDIR)/%.o: %.for
	@$(FC) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod
$(OBJDIR)/%.o: %.f90
	@$(FC) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

clean:
	@rm -f $(OBJDIR)/*.o

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
