# This is a makefile.

HOME = ./
PROG = capsid-souffle
OBJ = main.o initialize.o functions.o bead.o edge.o face.o subunit.o md.o analysis.o energies.o forces.o

# nanoHUB flags. 
nanoHUBCC = mpicxx -O3 -g -Wall -fopenmp -std=c++11 -Wunused-variable -Wunknown-pragmas
nanoHUBLFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_filesystem -lboost_system -lboost_program_options -lboost_mpi -lboost_serialization
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o

# BigRed2 flags. 
BigRed2CC = CC -O3 -g -Wall -fopenmp
BigRed2LFLAG = -L/N/soft/cle6/boost/gnu/1.72.0/lib -L/N/soft/cle6/gsl/2.6/lib -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
BigRed2CFLAG = -c
BigRed2OFLAG = -o

# General purpose flags.
CC = mpicxx -O3 -g -Wall -fopenmp -std=c++0x
LFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
CFLAG = -c
OFLAG = -o

all: $(PROG)

install:
	. /etc/environ.sh; use -e -r boost-1.62.0-mpich2-1.3-gnu-4.7.2; make CCF=nanoHUB all
	@echo "compiling the souffle code on Nanohub"
	cp -f $(PROG) $(BIN)

local-install: all
	@echo "creating output files folder: outfiles/"; mkdir $(HOME)outfiles

cluster-install:
	module swap PrgEnv-intel/6.0.5 PrgEnv-gnu && module load boost/gnu && module load gsl; make CCF=BigRed2 all
	@echo "Installing $(PROG) into $(HOME) directory on your computing cluster"; mkdir $(HOME)outfiles

$(PROG) : $(OBJ)
ifeq ($(CCF),BigRed2)	
	$(BigRed2CC) $(BigRed2OFLAG) $(PROG) $(OBJ) $(BigRed2LFLAG)
%.o : %.cpp
	$(BigRed2CC) -c $(BigRed2CFLAG) $< -o $@
else ifeq ($(CCF),nanoHUB)
	$(nanoHUBCC) $(nanoHUBOFLAG) $(PROG) $(OBJ) $(nanoHUBLFLAG)
%.o : %.cpp
	$(nanoHUBCC) -c $(nanoHUBCFLAG) $< -o $@
else
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)
%.o : %.cpp
	$(CC) -c $(CFLAG) $< -o $@	
endif
	
clean:
	rm -f *.o

dataclean: 
	rm -f outfiles/*.out outfiles/*.lammpstrj; rmdir $(HOME)outfiles
