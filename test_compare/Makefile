include ../GNUmakefile.chssi

#Compilers
FC=ftn
CC=cc
CPP=CC
SERIAL_CPP=CC

#Compiler flags
FFLAGS=-fastsse -Mcache_align -O3 -D_PORTGRP -DXT3 -D__fortran -D__fortran77 -DMPI2 -DNO_MPI_IO -mcmodel=medium -Mlarge_arrays
CFLAGS=-mcmodel=medium -Mlarge_arrays -fastsse -Mcache_align -O3 -D_PORTGRP -DC_SUFFIX -DCB_SUFFIX -DXT3 -DMPI2 -DNO_MPI_IO 
CPPFLAGS= -mcmodel=medium -Mlarge_arrays -Mcache_align -D_PORTGRP -DC_SUFFIX -DCB_SUFFIX -DXT3 -DMPI2 -DNO_MPI_IO 

#ACESLIBS variable will be replaced by the autoconf script with appropriate libraries
LIBS = -L.
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -lmpichf90

SIAL_COMPILER_LIBS=  -lsial -lsip_shared -laces2 -lrt -lpghpf2

LIB_DIRS:=-L../lib $(LIB_DIRS)

ACES2_LIB=../lib/libaces2.a
EXE=../bin/xtest_compare

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

MAIN_OBJS=$(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)

all: tabula_rasa $(EXE)

tabula_rasa: ;
	rm -f make.out make.err

$(EXE): xtest_compare
	cp xtest_compare $(EXE)
	rm -f make.out make.err

xtest_compare: $(ACES2_LIB) $(MAIN_OBJS)
	$(FC) $(FFLAGS) $(LFLAGS) -o xtest_compare \
	       $(MAIN_OBJS) \
	       $(LIB_DIRS) $(LIBS) && test -x xtest_compare  || rm -f xtest_compare

binclean: ;
	rm -f xtest_compare

ppclean clean distclean: % : tabula_rasa binclean
	rm -f *.o

rebuild relink: % : binclean $(EXE)

libclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

