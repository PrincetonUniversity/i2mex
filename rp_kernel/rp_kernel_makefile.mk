# Makefile to build librp_kernel.a and librp_kernel.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := rp_kernel
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

SRCF := $(wildcard *.for)
OBJF := $(SRCF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)

ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(FFLAGS) $(FPIC)


# Direct shared-library dependencies
L_RP_KERNEL_DEP := \
	-L$(LOCAL)/lib \
	-lrplot_mod \
	-lrplot_io \
	-lxdatmgr \
	-lcomput \
	-lportlib \
	-lsmlib \
	-L$(UFILES_ROOT)/lib \
	-lureadsub \
	-linterp_sub \
	-lmdstransp \
	-lr4smlib


# Runtime search paths
RPATH_RP_KERNEL := \
	-Wl,-rpath,'$$ORIGIN' \
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
	@$(FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) \
		$(L_RP_KERNEL_DEP) \
		$(RPATH_RP_KERNEL)
endif

include $(NAME)_depend.mk

$(OBJDIR)/%.o: %.for
	@$(FC) $(FPP) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

clean:
	@rm -f $(OBJDIR)/*.o

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
