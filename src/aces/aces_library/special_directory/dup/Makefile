include ../../../../../GNUmakefile.chssi

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

INCLUDE_DIRS := -I../../../../sia/include -I../aces2/include $(INCLUDE_DIRS)
LIB_DIRS:=-L../../../lib $(LIB_DIRS)

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

DUPLIB=../../../../../lib/libdup.a

all: tabula_rasa $(DUPLIB)

tabula_rasa: ;
	rm -f make.out make.err

$(DUPLIB): libdup.a
	cp libdup.a $(DUPLIB)
	rm -f make.out make.err

libdup.a: $(F_OBJS) $(f_OBJS) $(c_OBJS)
	ar $(ARFLAGS) libdup.a $(F_OBJS) $(f_OBJS) $(c_OBJS)

libclean: ;
	rm -f libdup.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(DUPLIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

