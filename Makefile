# This is a makefile.

HOME = ./

PROG = simulate_spinach_souffle

OBJ = main.o initialize.o functions.o bead.o edge.o face.o subunit.o md.o energies.o forces.o

CC = g++ -g -Wall -O3

LFLAG = -lgsl -lgslcblas

CFLAG = -c

OFLAG = -o

all: $(PROG)

install: all
	@echo "creating output files folder: outfiles/"; mkdir $(HOME)outfiles

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LIBS) $(LFLAG)

main.o:	md.h
md.o: initialize.h functions.h md.h energies.h forces.h
intialize.o: initialize.h rand_gen.h
function.o: functions.h bead.h subunit.h LJpair.h edge.h face.h
forces.o: forces.h bead.h LJpair.h subunit.h edge.h face.h functions.h
bead.o: bead.h edge.h
edge.o: edge.h bead.h face.h functions.h
energies.o: energies.h bead.h LJpair.h subunit.h edge.h face.h functions.h	
face.o: face.h bead.h edge.h functions.h
subunit.o: subunit.h bead.h

clean:
	rm -f *.o

dataclean: 
	rm -f outfiles/*.out outfiles/*.lammpstrj; rmdir $(HOME)outfiles
