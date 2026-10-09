# Makefile to build librplot_io.a and librplot_io.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := rplot_io
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := wkbuf_mod.mod
SRCM := $(MODS:%.mod=%.f90)
SF90 := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)
SRCF := $(SF90) $(SFOR)
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL) -I$(NETCDF_FORTRAN_HOME)/include
FFLG := $(FFLAGS)

L_RPLOT_IO_DEP := \
	-L$(LOCAL)/lib \
	-lrplot_mod \
	-lxdatmgr \
	-lportlib \
	-ltokyr \
	-lvaxonly \
	-lcomput \
	$(L_NETCDF) \
	$(L_MDSPLUS) \
	-L$(UFILES_ROOT)/lib \
	-lmdstransp \
	-lmds_sub \
	-lureadsub

RPATH_RPLOT_IO := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(UFILES_ROOT)/lib


.PHONY: all clean realclean clobber check depend

all: check $(NAME)

check:
	@mkdir -p $(LOCAL)/lib
	@mkdir -p $(OBJDIR)
	@mkdir -p $(LOCAL)/mod

$(NAME): $(OBJF)
	@echo "Building $(NAME) static library"
	@$(AR) $(LIBAR) $(OBJF)
ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@$(FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) \
		$(L_RPLOT_IO_DEP) \
		$(RPATH_RPLOT_IO)
endif

include $(NAME)_depend.mk

$(OBJDIR)/%.o: %.for
	@$(FC) $(FPP) $(FFLG) $< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mod

$(OBJDIR)/%.o: %.f90
	@$(FC) $(FPP) $(FFLG) $< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mod

clean:
	@rm -f $(OBJDIR)/*.o
	@cd $(LOCAL)/mod ; rm -f $(MODS)

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
