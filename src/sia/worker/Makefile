include ../../../GNUmakefile.chssi

#Compilers
FC=mpif77
CC=mpicc
CPP=mpicxx
SERIAL_CPP=g++

#Compiler flags
FFLAGS=-g
CFLAGS=
CPPFLAGS=

#ACESLIBS variable will be replaced by the autoconf script with appropriate libraries
LIBS = 
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= 

#MKLFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= 

SIAL_COMPILER_LIBS= 

INCLUDE_DIRS := -I../include -I../sip/include -I../../aces/aces_library/a2driver/aces2/include $(INCLUDE_DIRS)
LIB_DIRS:=-L../../../lib $(LIB_DIRS)

FSOURCE=$(wildcard *.F)
	F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

all rebuild relink: % : tabula_rasa $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS) tabula_rasa2

tabula_rasa: ;
	rm -f make.out make.err

tabula_rasa2: ;
	rm -f make.out make.err

ppclean clean distclean: % : tabula_rasa
	rm -f *.o

libclean binclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

