FC = gfortran
CC = gcc
FFLAGS = -O2 -Wall -ffixed-line-length-none
CFLAGS = -O2 -Wall -lm
LDLIBS = -L/usr/lib64 -lfftw3

FOBJS = whiten.o filter4.o sacio.c

whiten: $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o whiten $(LDLIBS)

clean:
	-rm *.o

