include ../GNUmakefile.chssi

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

A3ECP_LIB=../lib/liba3ecp.a

all: tabula_rasa $(A3ECP_LIB)

tabula_rasa: ;
	rm -f make.out make.err

$(A3ECP_LIB): liba3ecp.a
	cp liba3ecp.a $(A3ECP_LIB)

liba3ecp.a: $(F_OBJS) $(f_OBJS) $(c_OBJS)
	ar $(ARFLAGS) liba3ecp.a $(F_OBJS) $(f_OBJS) $(c_OBJS)

libclean: ;
	rm -f liba3ecp.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(A3ECP_LIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

