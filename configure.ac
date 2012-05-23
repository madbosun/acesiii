# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_HEADER([config.h])

# Checks for compilers
#If fortan compiler is not specified then try to autodetect it
if test x$FC = x; then
AC_PROG_FC(mpif77 mpif90 ftn)
fi
#If C compiler is not specified then try to autodetect it
if test x$CC = x; then
AC_PROG_CC(mpicc gcc cc)
fi
#If C++ compiler is not specified then try to autodetect it
if test x$CXX = x; then
AC_PROG_CXX(mpicxx g++ cc)
fi

AC_SUBST(SIAL_COMPILER_LIBS)
AC_SUBST(ACESLIBS)
AC_SUBST(ACESFLAGS)
AC_SUBST(MKLLIBS)
AC_SUBST(MKLFLAGS)
AC_SUBST(SERIAL_CPP)
AC_SUBST(CUDAPATH)

#This specifies the default C++ compiler in case no other compiler is found
SERIAL_CPP="g++"

#Settings for NVIDIA CUDA compiler and libraries
if test x$CUDAPATH = x; then
   AC_MSG_ERROR([Please specify the CUDA install path in variable CUDAPATH])
fi

#Settings for Intel MKL libraries on linux machine
AC_ARG_ENABLE([mkl],
    [AS_HELP_STRING([--enable-mkl], [Use Intel MKL Library])],
[
if test x$MKLLIBS = x; then
   AC_MSG_ERROR([Please specify all the required mkl libraries in variable MKLLIBS])
fi

if test x$MKLFLAGS = x; then
   AC_MSG_ERROR([Please specify mkl library location in variable  MKLFLAGS])
fi
#If intel mkl is used, use icpc as c++ compiler
SERIAL_CPP="icpc"
#Also set SIAL_COMPILER_LIBS accordingly
if test x$SIAL_COMPILER_LIBS = x; then
SIAL_COMPILER_LIBS="-lsial -lsip_shared -laces2 -lifcore"
fi
ACESLIBS=$MKLLIBS
ACESFLAGS=$MKLFLAGS
])


#Settings for GNU libraries on linux machine
AC_ARG_ENABLE([gnu],
    [AS_HELP_STRING([--enable-gnu], [Use GNU Library])],
[
if test x$GNULIBS = x; then
   AC_MSG_ERROR([Please specify all the required gnu libraries in variable GNULIBS])
fi

if test x$GNUFLAGS = x; then
   AC_MSG_ERROR([Please specify all the gnu library locations in variable  GNUFLAGS])
fi
#Set c++ compiler
SERIAL_CPP="g++"
#Also set SIAL_COMPILER_LIBS accordingly
if test x$SIAL_COMPILER_LIBS = x; then
SIAL_COMPILER_LIBS="-lsial -lsip_shared -laces2 -lgfortranbegin -lgfortran -ldl -lsip1 -lsip2"
fi
ACESLIBS=$GNULIBS
ACESFLAGS=$GNUFLAGS
])

#Settings for Cray XT3 system
AC_ARG_ENABLE([cray],
    [AS_HELP_STRING([--enable-cray], [Use ACESIII on Cray system])],
[
if test x$CRAYLIBS = x; then
   AC_MSG_ERROR([Please specify all the required libraries in variable CRAYLIBS])
fi

if test x$CRAYFLAGS = x; then
   AC_MSG_ERROR([Please specify all the library locations in variable CRAYFLAGS])
fi
#Set c++ compiler
SERIAL_CPP=$CXX
#Also set COMPILER_LIBS accordingly
if test x$SIAL_COMPILER_LIBS = x; then
SIAL_COMPILER_LIBS="-lsial -lsip_shared -lacescore -ltools -lsip1 -lsip2"
fi
ACESLIBS=$CRAYLIBS
ACESFLAGS=$CRAYFLAGS
])


# Checks for libraries.
# No checks yet

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/main/Makefile])
AC_CONFIG_FILES([src/sial_compiler/compiler/Makefile])
AC_CONFIG_FILES([src/sial_compiler/sial/Makefile])
AC_CONFIG_FILES([src/aces/aces_sial/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/special_directory/oed/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/special_directory/dup/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/special_directory/nlo/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/special_directory/erd/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/special_directory/direct_integral/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/a2driver/symcor/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/a2driver/aces2/Makefile])
AC_CONFIG_FILES([src/aces/aces_library/a2driver/geopt/Makefile])
AC_CONFIG_FILES([src/sia/framelib/Makefile])
AC_CONFIG_FILES([src/sia/manager/Makefile])
AC_CONFIG_FILES([src/sia/sip/aces_instructions/Makefile])
AC_CONFIG_FILES([src/sia/sip/sip_instructions/Makefile])
AC_CONFIG_FILES([src/sia/worker/Makefile])
AC_CONFIG_FILES([src/sia/sip_shared/Makefile])
AC_CONFIG_FILES([test_compare/Makefile])

AC_OUTPUT
