# This is a makefile.

HOME = ./
BIN = ../bin
PROG = capsid-souffle
OBJ = main.o initialize.o functions.o bead.o edge.o face.o subunit.o md.o analysis.o energies.o forces.o

# nanoHUB flags. 
nanoHUBCC = mpicxx -O3 -g -Wall -fopenmp -std=c++11 -Wunused-variable -Wunknown-pragmas
nanoHUBLFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_filesystem -lboost_system -lboost_program_options -lboost_mpi -lboost_serialization
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o

# BigRed2 flags. 
BigRed2CC = CC -O3 -g -Wall -fopenmp -std=c++11 -Wunused-variable -Wunknown-pragmas
BigRed2LFLAG = -L/N/soft/cle7/boost/gnu/1.76.0/lib -L/N/soft/cle7/gsl/gnu/2.7/lib -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -lboost_program_options -lboost_mpi -lboost_serialization
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
	module swap PrgEnv-intel PrgEnv-gnu && module load boost/gnu && module load gsl/gnu/2.7; make CCF=BigRed2 all
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
	rm -f $(PROG)

dataclean: 
	rm -f $(BIN)/outfiles/*.out $(BIN)/outfiles/*.lammpstrj
	rmdir $(BIN)/outfiles
	
distclean: clean
	rm -f $(BIN)/$(PROG)
