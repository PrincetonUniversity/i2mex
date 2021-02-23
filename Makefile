#-*- makefile -*-
#GNU Makefile to build Library libi2mex.a 
#
#12/04/00 C. Ludescher
#Modified 01/31/20 J. Breslau
#module load intel; module load ntcc/stable lapack mdsplus

OBJ=.
FC90=ifort
CC = icc

#FFLAGS += -check underflow -check overflow -check bounds \
#-warn argument_checking

FFLAGS += -O

.PHONY: all clean

ARC =  $(OBJ)/lib/libi2mex.a
OBJDIR = $(OBJ)/obj

# f90 modules
Mnams = cont_mod.mod i2mex_mod.mod freeqbe_mod.mod
Mobjs = cont_mod.o i2mex_mod.o freeqbe_mod.o
MODS = $(foreach file,$(Mobjs),$(OBJDIR)/$(file))

# library members
ALLMEM = $(subst .f90,.o, $(wildcard *.f90))
DFMEM =$(filter-out mex2eqs.o, $(filter-out drive.o, $(ALLMEM)))
FMEM =$(filter-out cont_mod.o, $(filter-out freeqbe_mod.o, $(DFMEM)))
MEM = $(foreach m,$(FMEM),$(OBJDIR)/$(m))

ifdef NO_TRXPLIB
      TRXPLIB= -lgeneric_dummy
else
      TRXPLIB=-ltrxplib -ltrread -ltr_getnl -lrp_kernel -lrplot_io \
	-lmds_sub -lmdstransp -lxdatmgr -linterp_sub
endif

LDLIBS2 = -L$(OBJ)/lib -L$(NTCC_HOME)/lib -llsode -llsode_linpack $(TRXPLIB) \
 -lold_xplasma -lxplasma_debug -lxplasma2 -lgeqdsk_mds -lr8bloat \
 -lsmlib -lmdstransp -L$(MDSPLUS)/lib -lMdsLib -lnscrunch -lfluxav \
 -ltrgraf  -L$(PSPLINE_HOME)/lib -lpspline \
 -L${NETCDF_FORTRAN_HOME}/lib -lnetcdf -lnetcdff \
 -L${LAPACK_HOME}/lib -llapack -lrefblas -lezcdf \
 -lmclib -lureadsub -lcomput -lvaxonly -llsode -llsode_linpack \
 -lelvislib -lsg -ljc -lportlib 

LDLIBS = -L$(OBJ)/lib -L$(NTCC_HOME)/lib -llsode -llsode_linpack $(TRXPLIB) \
 -lold_xplasma -lxplasma2 -lgeqdsk_mds -lr8bloat -lmdstransp \
 -L$(MDSPLUS)/lib -lMdsLib -lnscrunch -lsmlib -lfluxav -L$(PSPLINE_HOME)/lib -lpspline \
 -L${NETCDF_FORTRAN_HOME}/lib -lnetcdf -lnetcdff \
 -L${LAPACK_HOME}/lib -llapack -lrefblas -lezcdf -lmclib \
 -lcomput -lvaxonly -lportlib

libs: chkdirs $(ARC)

all: libs exec
	@echo done 

$(OLDLIB): timestamp.inf
	@echo "--- DETECTED i2mex.a source update"
	@echo "--- Re-Making $(ARC)"
	@$(MAKE) libs

chkdirs:
	@test -d $(OBJ)/lib || mkdir -p $(OBJ)/lib
	@test -d $(OBJ)/mod || mkdir -p $(OBJ)/mod
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)
	@test -d $(OBJ)/bin || mkdir -p $(OBJ)/bin

# compile f90
$(OBJDIR)/%.o: %.f90
	$(FC90) $(FFLAGS) $(MODFLAGS) -I./ -I$(NTCC_HOME)/mod -I$(PSPLINE_HOME)/include -c -o $(OBJDIR)/$*.o $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

$(ARC): $(MODS) $(MEM)
	ar -r $(ARC) $(MODS) $(MEM)
	cp -u *.mod $(OBJ)/mod

#JB exec: $(OBJ)/test/i2mex $(OBJ)/test/mex2eqs makelink 
exec: $(OBJ)/bin/i2mex $(OBJ)/bin/mex2eqs

$(OBJ)/bin/i2mex: $(OBJDIR)/drive.o $(ARC)
	$(FC90) $(LDFLAGS) -o $@ $< $(ARC) $(LDLIBS)

$(OBJ)/bin/mex2eqs: $(OBJDIR)/mex2eqs.o readline.o $(ARC)
	$(FC90) $(LDFLAGS) -o $@ $< readline.o $(ARC) $(LDLIBS2)

clean:
	@rm -f $(OBJDIR)/* readline.o
	@rm -f $(OBJ)/mod/*
	@rm -f *.mod
	@rm -f $(OBJ)/bin/*
	@rm -f $(ARC)
