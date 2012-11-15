include ../../../../GNUmakefile.chssi

.SUFFIXES: .cu

#Compilers
FC=mpif77
CC=mpicc
CPP=mpicxx
NVCC=nvcc
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

SIAL_COMPILER_LIBS= 

INCLUDE_DIRS := -I../../include -I../include -I../../../aces/aces_library/a2driver/aces2/include $(INCLUDE_DIRS)
LIB_DIRS:= -L../../../../lib $(LIB_DIRS)

SIP_LIB=../../../../lib/libsip2.a

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)
CU_SOURCE=$(wildcard *.cu)
CU_OBJS:=$(CU_SOURCE:.cu=.o)


all: tabula_rasa $(SIP_LIB)

tabula_rasa: ;
	rm -f make.out make.err

$(SIP_LIB): libsip2.a
	cp libsip2.a $(SIP_LIB)
	rm -f make.out make.err

libsip2.a: $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS) $(CU_OBJS)
	ar $(ARFLAGS) libsip2.a $(F_OBJS) $(f_OBJS) $(c_OBJS) \
	   $(CPP_OBJS) $(CU_OBJS)

libclean: ;
	rm -f libsip2.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(SIP_LIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cu.o:
	$(NVCC) -arch=sm_20 -Xcompiler $(CPPFLAGS) -c $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err
