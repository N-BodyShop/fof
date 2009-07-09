#
# Makefile for fof.
#
# JPG 5/29/09: This has been modified to include OpenMP
#              capability.  The algorithm could be improved,
#              but at least it makes some use of extra cores
#              that might be lying around.  If not compiled
#              with OpenMP, the code should be absolutely
#              identical to the original serial version.
#              If compiled with OpenMP and run on 1 core,
#              there will be a slight performance hit due
#              to overhead from the thread locking infrastructure.


#CFLAGS	=	-O2
# Compile for Linux gcc
#CFLAGS  =   -O3 -funroll-loops
#CFLAGS = -O3
# Compile Intel statically linked
#CC = icc
#CFLAGS = -O3 -ipo -static
# Compile Intel Openmp
CC = icc
CFLAGS = -O3 -ipo -openmp -openmp-report2


LIBS	=   -lm

default:	fof clean

clean:
	rm -f *.o

fof: main.o kd.o
	$(CC) $(CFLAGS) -o fof main.o kd.o $(LIBS)

main.o: main.c kd.h
	$(CC) $(CFLAGS) -c main.c

kd.o: kd.c kd.h tipsydefs.h
	$(CC) $(CFLAGS) -c kd.c



