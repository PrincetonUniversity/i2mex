# Makefile to build libtr_getnl.a and libtr_getnl.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := tr_getnl
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

MODS := tr_getnl.mod
SRCM := $(MODS:%.mod=%.f90)
SRCF := $(SRCM) $(filter-out $(SRCM), $(wildcard *.f90))
OBJF := $(SRCF:%.f90=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)
FFLG := $(FFLAGS)

L_TR_GETNL_DEP := \
	-L$(LOCAL)/lib \
	-lportlib \
	-lrplot_io

RPATH_TR_GETNL := \
	-Wl,-rpath,'$$ORIGIN'


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
		$(L_TR_GETNL_DEP) \
		$(RPATH_TR_GETNL)
endif

include $(NAME)_depend.mk

$(OBJDIR)/%.o: %.f90
	@$(FC) $(FFLG) $< \
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
