#   
#   Main Makefile for aces  application.
#

export
CHSSI_EXE=xchssi

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
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -lecp -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -lmpichf90

SIAL_COMPILER_LIBS=  -lsial -lsip_shared -laces2 -lrt -lpghpf2


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


