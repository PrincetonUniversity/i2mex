# Makefile to build libtrxplib.a
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := trxplib
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := trx_bxtr_options.mod trx_module.mod trxplib_ps_options.mod
SRCM := $(MODS:%.mod=%.f90)
SF90 := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)
SRCF := $(SF90) $(SFOR)
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL) -I$(PLASMA_STATE_ROOT)/mod

L_TRXPLIB_DEP := \
	-L$(LOCAL)/lib \
	-lxdatmgr \
	-lrp_kernel \
	-ltrread \
	-lrplot_io \
	-lsplitn \
	-ltr_getnl \
	-lechmod_iolib \
	-ltrdatbuf_lib \
	-L$(PLASMA_STATE_ROOT)/lib \
	-lplasma_state \
	-lold_xplasma \
	-lxplasma2 \
	-lps_xplasma2

RPATH_TRXPLIB := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(PLASMA_STATE_ROOT)/lib


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
	@$(FC) -shared $(LDFLAGS) -o $(LIBSO) $(OBJF) \
	$(L_TRXPLIB_DEP) \
	$(RPATH_TRXPLIB)
#-L$(PLASMA_STATE_ROOT)/lib -lplasma_state \
#	-L$(LOCAL)/lib -ltrdatbuf_lib 
endif

include $(NAME)_depend.mk
$(OBJDIR)/%.o: %.for
	@$(FC) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod
$(OBJDIR)/%.o: %.f90
	@$(FC) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod ; rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
