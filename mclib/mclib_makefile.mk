# Makefile to build libmclib.a and libmclib.so
#

# Include system dependent flags
include $(COMPILER_FLAGS)

NAME := mclib
OBJDIR := $(LOCAL)/obj/$(NAME)
LIBAR := $(LOCAL)/lib/lib$(NAME).a
LIBSO := $(LOCAL)/lib/lib$(NAME).so

SRCF := $(wildcard *.for)
OBJF := $(SRCF:%.for=$(OBJDIR)/%.o)

FDEF := $(FDEFS)
FINC := $(FINCL)
FFLG := $(FFLAGS)

L_MCLIB_DEP := \
	-L$(LOCAL)/lib

RPATH_MCLIB := \
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
		$(L_MCLIB_DEP) \
		$(RPATH_MCLIB)
endif

include $(NAME)_depend.mk

$(OBJDIR)/%.o: %.for
	@$(FC) $(FPP) $(FFLG) $< \
		-o $@ \
		$(FDEF) \
		$(FINC) \
		$(MFLAG) $(LOCAL)/mod

clean:
	@rm -f $(OBJDIR)/*.o

realclean: clean
	@rm -f $(LIBAR) $(LIBSO)

clobber: realclean

depend:
	@makedepf90 -b OBJDIR $(SRCF) > $(NAME)_depend.mk
	@sed -i 's/OBJDIR/\$$\(OBJDIR\)/g' $(NAME)_depend.mk
