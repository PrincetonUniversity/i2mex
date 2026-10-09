# Makefile to build libtrdatbuf_lib.a
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := trdatbuf_lib
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := trdatbuf_obj.mod trdatbuf_aux.mod tdb_static.mod tdbsub_uts.mod \
	trdatbuf_iface.mod trdatbuf_intmod.mod trdatbuf_module.mod
SRCM := $(MODS:%.mod=%.f90)
SF90 := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
SFOR := $(wildcard *.for)
SRCF := $(SF90) $(SFOR)
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)
OBJF := $(OBJF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL) -I$(NETCDF_FORTRAN_HOME)/include -I$(UFILES_ROOT)/mod

ifeq ($(MAKE_SO),1)
	FPIC := -fPIC
else
	FPIC :=
endif

FFLG := $(FFLAGS) $(FPIC)

L_TRDATBUF_LIB_DEP := \
	-L$(LOCAL)/lib \
	-lportlib \
	-lcomput \
	$(L_PSPLINE) \
	$(L_UFILES) \
	$(L_NETCDF)

RPATH_TRDATBUF_LIB := \
	-Wl,-rpath,'$$ORIGIN' \
	-Wl,-rpath,$(PSPLINE_HOME)/lib \
	-Wl,-rpath,$(UFILES_ROOT)/lib \
	-Wl,-rpath,$(NETCDF_FORTRAN_HOME)/lib \
	-Wl,-rpath,$(NETCDF_C_HOME)/lib64 \
	-Wl,-rpath,$(MDSPLUS_ROOT)/lib


.PHONY: clean realclean clobber check

all: check $(NAME)

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
		$(L_TRDATBUF_LIB_DEP) \
		$(RPATH_TRDATBUF_LIB)
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
