# Makefile to build libvaxonly.a
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := vaxonly
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

SRCF := $(wildcard *.for)
OBJF := $(SRCF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)
FFLG := $(FFLAGS)

L_VAXONLY_DEP := \
	-L$(LOCAL)/lib \
	-lportlib \
	$(L_NETCDF)

RPATH_VAXONLY := \
	-Wl,-rpath,'$$ORIGIN'


.PHONY: all vaxonly clean realclean clobber check depend

all: check vaxonly

vaxonly: check $(OBJF)
	@echo "Building $(NAME) static library"
	@$(AR) $(LIBAR) $(OBJF)
ifeq ($(MAKE_SO),1)
	@echo "Building $(NAME) shared library"
	@$(FC) -shared $(LDFLAGS) \
		-o $(LIBSO) \
		$(OBJF) \
		$(L_VAXONLY_DEP) \
		$(RPATH_VAXONLY)
endif

check:
	@test -d $(LOCAL)/lib || mkdir -p $(LOCAL)/lib
	@test -d $(LOCAL)/obj/$(NAME) || mkdir -p $(LOCAL)/obj/$(NAME)
	@test -d $(LOCAL)/mod || mkdir -p $(LOCAL)/mod

include $(NAME)_depend.mk

$(OBJDIR)/%.o: %.for
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
