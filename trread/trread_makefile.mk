# Makefile to build libtrread.a and libtrread.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := trread
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := tconnect_mod.mod
SRCM := $(MODS:%.mod=%.f90)
SF90 := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)
SRCF := $(SF90) $(SFOR)
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

SRCC := $(wildcard *.c)
OBJC := $(SRCC:%.c=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)

CDEF := $(CDEFS)
CINC := $(CINCL) -I../incl_cpp

ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(FFLAGS) $(FPIC)
CFLG := $(CFLAGS) $(FPIC)


# Direct shared-library dependencies
L_TRREAD_DEP := \
	-L$(LOCAL)/lib \
	-lrplot_mod \
	-lrp_kernel \
	-lrplot_io \
	-ltr_getnl \
	-lxdatmgr \
	-lcomput \
	-lvaxonly \
	-lportlib \
	-L$(UFILES_ROOT)/lib \
	-lmdstransp \
	-linterp_sub


# Runtime search paths
RPATH_TRREAD := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(UFILES_ROOT)/lib


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
	@$(FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) $(OBJC) \
		$(L_TRREAD_DEP) \
		$(RPATH_TRREAD)
endif

include $(NAME)_depend.mk

$(OBJDIR)/%.o: %.for
	@$(FC) $(FPP) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

$(OBJDIR)/%.o: %.f90
	@$(FC) $(FFLG) $< -o $@ $(FDEF) $(FINC) $(MFLAG) $(LOCAL)/mod

$(OBJDIR)/%.o: %.c
	@$(CC) $(CFLG) $< -o $@ $(CDEF) $(CINC)

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod ; rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
