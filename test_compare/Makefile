include ../GNUmakefile.chssi

#Compilers
FC=mpif77
CC=mpicc
CPP=mpicxx
SERIAL_CPP=icpc

#Compiler flags
FFLAGS=-D__fortran -D__fortran77 -DMPIF2C -Zp8 -zero -traceback -mcmodel=medium -DMPI2 -shared-intel -O2
CFLAGS=-DMPIF2C -DC_SUFFIX -mcmodel=medium -DMPI2 -shared-intel -O2  -DCB_SUFFIX
CPPFLAGS=-DMPIF2C -DC_SUFFIX -DCB_SUFFIX -mcmodel=medium -DMPI2 -shared-intel -O2

#ACESLIBS variable will be replaced by the autoconf script with appropriate libraries
LIBS = -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -Wl,--end-group
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -L.

SIAL_COMPILER_LIBS= -lsial -lsip_shared -laces2  -lifcore

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

