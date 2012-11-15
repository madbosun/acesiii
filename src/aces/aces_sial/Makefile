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

INCLUDE_DIRS := -I../../sia/include -I../../sia/sip/include -I../aces_library/a2driver/aces2/include $(INCLUDE_DIRS)
SIAL_COMPILER_LIBS= 


SIAL_COMPILER=../../../bin/sial
LIB_SIALINST=../../../lib/libsialinst.a
SIO_DIR=../../../bin/sio

#Search for all the .sial files in current directory as well as in sub directories upto 3 levels deep
SIAL_SOURCE := $(wildcard *.sial)
SIAL_SOURCE += $(wildcard */*.sial)
SIAL_SOURCE += $(wildcard */*/*.sial)
SIAL_SOURCE += $(wildcard */*/*/*.sial)
SIO_OBJS=$(SIAL_SOURCE:.sial=.sio)
#Search for all the fortran files in current directory as well as in sub directories upto 3 levels deep
FSOURCE :=$(wildcard *.F)
FSOURCE +=$(wildcard */*.F)
FSOURCE +=$(wildcard */*/*.F)
FSOURCE +=$(wildcard */*/*/*.F)
F_OBJS =$(FSOURCE:.F=.o)
fSOURCE :=$(wildcard *.f)
fSOURCE +=$(wildcard */*.f)
fSOURCE +=$(wildcard */*/*.f)
fSOURCE +=$(wildcard */*/*/*.f)
f_OBJS :=$(fSOURCE:.f=.o)
#Search for all the C/C++ files in current directory as well as in sub directories upto 3 levels deep
cSOURCE :=$(wildcard *.c)
cSOURCE +=$(wildcard */*.c)
cSOURCE +=$(wildcard */*/*.c)
cSOURCE +=$(wildcard */*/*/*.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE :=$(wildcard *.cpp)
CPP_SOURCE +=$(wildcard */*.cpp)
CPP_SOURCE +=$(wildcard */*/*.cpp)
CPP_SOURCE +=$(wildcard */*/*/*.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

LD_LIBRARY_PATH :=/opt/intel/mkl/11.1.064/lib/em64t:${LD_LIBRARY_PATH}
all: tabula_rasa $(LIB_SIALINST) $(SIO_OBJS) $(SIO_DIR)/sial_config $(SIO_DIR)/default_jobflows

tabula_rasa: ;
	rm -f make.out make.err

binclean clean distclean: % : ;
	find . -type f -name "*.sio" -exec rm -f '{}' \;
	find . -type f -name "*.o" -exec rm -f '{}' \;
rebuild relink: % : binclean all

libclean ppclean archive: % : ;

$(LIB_SIALINST): libsialinst.a
	cp libsialinst.a $(LIB_SIALINST)
	rm -f make.out make.err

libsialinst.a: $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)
	ar $(ARFLAGS) libsialinst.a $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)

$(SIO_DIR)/sial_config: sial_config
	cp sial_config $(SIO_DIR)

$(SIO_DIR)/default_jobflows: default_jobflows
	cp default_jobflows $(SIO_DIR)

.SUFFIXES: .sial .sio

.sial.sio:
	$(SIAL_COMPILER) $< && cp $@ $(SIO_DIR)

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< -o $@  1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< -o $@ 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< -o $@  1>>make.out 2>>make.err

.cpp.o:
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< -o $@ 1>>make.out 2>>make.err

