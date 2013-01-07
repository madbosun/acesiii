#   
#   Main Makefile for aces  application.
#

export
CHSSI_EXE=xchssi

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
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -lecp -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -L.

SIAL_COMPILER_LIBS= -lsial -lsip_shared -laces2  -lifcore


FILES:= src src/aces/aces_sial test_compare tests
TARGET_DIRS:=$(shell for dir in $(FILES); \
                     do test -f $$dir/Makefile && echo $$dir; \
                     done)

all binclean libclean ppclean clean distclean: % : ;
	@for dir in $(TARGET_DIRS) ; \
	 do $(MAKE) -C $$dir $@ || exit 1 ; \
	 done

relink: binclean all

rebuild: libclean all

archive:


