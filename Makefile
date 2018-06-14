# This is a makefile.

PROG = simulate_spinach_souffle

OBJ = main.o initialize.o functions.o bead.o edge.o face.o unit.o md.o

CC = g++ -g -Wall

LFLAG = -lgsl -lgslcblas

CFLAG = -c

OFLAG = -o

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LIBS) $(LFLAG)

main.o:	md.h
md.o: initialize.h functions.h md.h
intialize.o: initialize.h rand_gen.h
function.o: functions.h bead.h unit.h LJpair.h edge.h face.h
bead.o: bead.h edge.h
edge.o: edge.h bead.h face.h functions.h
face.o: face.h bead.h edge.h functions.h
unit.o: unit.h bead.h

clean:
	rm -f *.o
