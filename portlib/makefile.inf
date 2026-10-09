#GNU Makefile to build portlib.a 
#        and test driver porttest
#
#01/08/09 C. Ludescher-Furth
#

libs:
	$(MAKE) -f Make2 libs
ifdef MPI_ENABLE
	$(MAKE) -f Make2 MPI_ACTIVATE=1 libs
endif

all: libs exec

exec:
	$(MAKE) -f Make2 exec
ifdef MPI_ENABLE
	$(MAKE) -f Make2 MPI_ACTIVATE=1 exec
endif

clean:
	$(MAKE) -f Make2 clean
ifdef MPI_ENABLE
	$(MAKE) -f Make2 MPI_ACTIVATE=1 clean
endif

realclean:
	$(MAKE) -f Make2 realclean
ifdef MPI_ENABLE
	$(MAKE) -f Make2 MPI_ACTIVATE=1 realclean
endif

install:
	$(MAKE) -f Make2 install
ifdef MPI_ENABLE
	$(MAKE) -f Make2 MPI_ACTIVATE=1 install
endif

uninstall:
	$(MAKE) -f Make2 uninstall
ifdef MPI_ENABLE
	$(MAKE) -f Make2 MPI_ACTIVATE=1 uninstall
endif

