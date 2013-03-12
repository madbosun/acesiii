include ../../../GNUmakefile.chssi

#Compilers
FC=mpif77
CC=mpicc
CPP=mpicxx
SERIAL_CPP=icpc

#Compiler flags
FFLAGS=-mcmodel=medium -D__fortran -D__fortran77 -DMPIF2C -Zp8 -zero -traceback -DMPI2 -shared-intel -O2
CFLAGS=-mcmodel=medium -DMPIF2C -DC_SUFFIX  -DMPI2 -shared-intel -O2  -DCB_SUFFIX
CPPFLAGS=-mcmodel=medium -DMPIF2C -DC_SUFFIX -DCB_SUFFIX  -DMPI2 -shared-intel -O2

#ACESLIBS variable will be replaced by the autoconf script with appropriate libraries
LIBS = -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -Wl,--end-group
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -L.

SIAL_COMPILER_LIBS= -lsial -lsip_shared -laces2  -lifcore

INCLUDE_DIRS := -I../../sia/include -I../sial -I../../aces/aces_library/a2driver/aces2/include/ $(INCLUDE_DIRS)
LIB_DIRS:=-L../../../lib $(LIB_DIRS)

SIAL_INCLUDES=$(INCLUDE_DIRS) -I../../aces/aces_sial/
SIAL_LIB=../../../lib/libsial.a
SIAL_COMPILER=../../../bin/sial
SIP_SHARED_LIB=../../../lib/libsip_shared.a

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

all: tabula_rasa $(SIAL_COMPILER) 

tabula_rasa: ;
	rm -f make.out make.err

$(SIAL_COMPILER): sial
	cp sial $(SIAL_COMPILER)
	rm -f make.out make.err
sial: $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS) $(SIAL_LIB) $(SIP_SHARED_LIB) 
	$(SERIAL_CPP) $(CPPFLAGS) $(SIAL_INCLUDES) $(LIB_DIRS) -o sial \
	$(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS) $(SIAL_COMPILER_LFLAGS) \
	$(SIAL_COMPILER_LIBS) \
	&& test -x sial || rm -f sial

binclean: ;
	rm -f sial

ppclean clean distclean: % : tabula_rasa binclean
	rm -f *.o

rebuild relink: % : binclean $(SIAL_COMPILER)

libclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

