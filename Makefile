#-*- makefile -*-
#GNU Makefile to build Library libi2mex.a 
#
#12/04/00 C. Ludescher

#JB ifneq ($(MAKELEVEL),0)
#JB # if OBJ was defined for main Makefile one level down
#JB ifeq ("${OBJ}",".")
#JB 	override OBJ=..
#JB endif
#JB endif

OBJ=.
FC90=ifort

# define system dependent flags, etc.
#JB -include  ../share/Make.local
#JB include ../share/Make.flags



ifdef DEBUG
        FFLAGS := $(DFFLAGS)
endif

#FFLAGS += -check underflow -check overflow -check bounds \
#-warn argument_checking

.PHONY: all install clean realclean

ARC =  $(OBJ)/lib/libi2mex.a
OBJDIR = $(OBJ)/obj

# don't rebuild library if using installed one in $PREFIX 
ifeq ($(MAKELEVEL),0)
	THISLIB=$(ARC)
endif

srcdir = $(shell pwd)
DATA   = $(shell ls *.cdf INP1* *.m *.net)

# f90 modules
Mnams = cont_mod.mod i2mex_mod.mod freeqbe_mod.mod
Mobjs = cont_mod.o i2mex_mod.o freeqbe_mod.o
#JB MODS = $(foreach file,$(Mobjs),$(ARC)($(file)))
MODS = $(foreach file,$(Mobjs),$(OBJDIR)/$(file))

ifeq ($(MODUP),Y)
 MODS0=$(foreach m,$(Mnams),$(shell  echo $(m) | tr 'a-z' 'A-Z'))
 MODULES=$(foreach m,$(MODS0),$(subst .MOD,.$(MODEXT),$(m)))
else
 MODULES = $(foreach m,$(Mnams),$(subst .mod,.$(MODEXT),$(m)))
endif 

# library members
ALLMEM = $(subst .f90,.o, $(wildcard *.f90))
#JB FMEM =$(filter-out mex2eqs.o, $(filter-out drive.o, $(ALLMEM)))
FMEM =$(filter-out cont_mod.o, $(filter-out freeqbe_mod.o, $(ALLMEM)))
#JB MEM = $(foreach m,$(FMEM),$(ARC)($(m)))
MEM = $(foreach m,$(FMEM),$(OBJDIR)/$(m))

ifdef NO_TRXPLIB
      TRXPLIB= -lgeneric_dummy
#      MDSLIB= -lmds_dummy
#      MDSTR= -lmdstransp
else
      TRXPLIB=-ltrxplib -ltrread -ltr_getnl -lrp_kernel -lrplot_io \
	-lmds_sub -lmdstransp -lxdatmgr -linterp_sub
#  ifdef MDSLIB
#      MDSLIB:= $(MDSLIB) -lTreeShr -lMdsShr -lTdiShr
#      MDSTR=-lmdstransp 
#  else
#      MDSLIB = -lmds_dummy
#      MDSTR=
#  endif
endif

LDLIBS2 = -L$(OBJ)/lib $(LLOC) -lesc -llsode -llsode_linpack $(TRXPLIB) \
	  -lold_xplasma -lxplasma_debug -lxplasma2 -lgeqdsk_mds -lr8bloat \
	  -lsmlib -lmdstransp -lnscrunch -lfluxav \
	  -ltrgraf -lpspline -lezcdf \
	  -lmclib -lureadsub -lcomput -lvaxonly -llsode -llsode_linpack \
	  -lelvislib -lsg -ljc -lportlib 

LDLIBS = -L$(OBJ)/lib $(LLOC) -lesc -llsode -llsode_linpack $(TRXPLIB) \
         -lold_xplasma -lxplasma2 -lgeqdsk_mds -lr8bloat -lmdstransp\
	 -lnscrunch -lsmlib -lfluxav -lpspline -lezcdf -lmclib \
	-lcomput -lvaxonly -lportlib  

libs: chkdirs $(ARC)

#JB all: libs exec
all: libs
	@echo done 

$(OLDLIB): timestamp.inf
	@echo "--- DETECTED i2mex.a source update"
	@echo "--- Re-Making $(ARC)"
	@$(MAKE) libs

chkdirs:
	@test -d $(OBJ)/lib || mkdir -p $(OBJ)/lib
	@test -d $(OBJ)/mod || mkdir -p $(OBJ)/mod
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)
#JB 	@test -d $(OBJ)/test || mkdir $(OBJ)/test

# compile f90
$(OBJDIR)/%.o: %.f90
	$(FC90) $(FFLAGS) $(MODFLAGS) -I./ -I$(NTCC_HOME)/mod -I$(PSPLINE_HOME)/include -c -o $(OBJDIR)/$*.o $<
	ar -r $(ARC) $(OBJDIR)/$*.o

$(ARC): $(MODS) $(MEM)
	cp -u *.mod $(OBJ)/mod
#JB 	$(RANLIB) $@ > /dev/null

#JB exec: $(OBJ)/test/i2mex $(OBJ)/test/mex2eqs makelink 

#JB makelink: 
#JB 	@for i in $(DATA); do \
#JB 	(cd $(OBJ)/test; test -f $$i || ln -s $(srcdir)/$$i $$i;) done 



install: installdirs installlib

installlib: 
	umask 133; cp $(OBJ)/lib/libi2mex.a $(LIBDIR)
	$(foreach file,$(MODULES),$(shell cp $(MDIR)/$(file) $(MODDIR)/))
#	umask 133; cp get_c.3 $(MANDIR)/man3

installdirs:
	@test -d $(MODDIR) ||  mkdir -p $(MODDIR)
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
#	@test -d $(MANDIR)/man3 || mkdir -p $(MANDIR)/man3


clean:
	@rm -f $(OBJDIR)/*
	@rm -f $(OBJ)/mod/*
#JB 	@(cd $(OBJ)/test; rm -f i2mex* )
	@if test -d $(MDIR); then \
	  (cd $(MDIR); rm -f $(MODULES)); fi

realclean: clean
	@rm -f $(ARC)
#JB 	@(cd $(OBJ)/test; rm -f $(DATA))

uninstall:
	$(foreach m,$(MODULES),$(shell rm -f $(MODDIR)/$(m) ))
	rm -f $(LIBDIR)/libi2mex.a




