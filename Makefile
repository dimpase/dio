#########################################################################
##
## Makefile for Hilbert basis package, version 0.1 (experimental)
## Copyright (C) 1995	Dmitrii V. Pasechnik, 
##			RIACA, Amsterdam, The Netherlands
##
## to compile all the programs, type 'make'
## it will create files oneeq.exe and systeq.exe, respectively for
## solving one (maybe inhomogeneous) linear Diophantine equation,
## and a system of homogeneous linear Diophantine equations.
##
## for the time being, you cannot use the standard Unix C-compiler (cc)
## 
## you can also create only one of those programs by typing 'make oneeq'
## (or 'make systeq', respectively)
##
## typing 'make clean' deletes the object files, whereas
## 'make allclean' deletes also *.exe files.
##
CC=gcc
CFLAGS=-O -g # -pg
EXT=.exe

all : oneeq systeq

oneeq : oneeq.o hb.o hb_uts.o hb.h
	$(CC) $(CFLAGS) -o oneeq$(EXT) oneeq.o hb.o hb_uts.o

systeq : systeq.o hb.o hb_uts.o hbs.o hb.h
	$(CC) $(CFLAGS) -o systeq$(EXT) systeq.o hb.o hb_uts.o hbs.o 

oneeq.o : oneeq.c hb.h
	$(CC) $(CFLAGS) -c oneeq.c

systeq.o : systeq.c hb.h
	$(CC) $(CFLAGS) -c systeq.c

hb.o : hb.c hb.h
	$(CC) $(CFLAGS) -c hb.c

hbs.o : hbs.c hb.h
	$(CC) $(CFLAGS) -c hbs.c

hb_uts.o : hb_uts.c hb.h
	$(CC) $(CFLAGS) -c hb_uts.c

clean :
	rm oneeq.o systeq.o hb.o hbs.o hb_uts.o 

allclean :
	rm oneeq.o systeq.o hb.o hbs.o hb_uts.o oneeq$(EXT) systeq$(EXT)

