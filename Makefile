#   
#   Main Makefile for aces  application.
#

export
CHSSI_EXE=xchssi

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

SIAL_COMPILER_LIBS= 


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


