#
# Makefile for fof.
#
CFLAGS	=	-O2
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



