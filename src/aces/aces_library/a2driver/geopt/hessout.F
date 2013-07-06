      SUBROUTINE HESSOUT(X,IR,IC,NSYM)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "coord.com"
#include "flags.h"
#include "jodaflags.com"
C
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
#include "cbchar.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
#include "io_units.par"
      DIMENSION X(IR,IC)
C
      JX=1
 1    ITOP=MIN(JX+4,IR)
      
      IF(NSYM.EQ.1)
     $  WRITE(LuOut,90)(VARNAM(ISQUASH(NOPTI(IX))),IX=JX,ITOP)
      IF(NSYM.EQ.0)WRITE(LuOut,90)(VARNAM(ISQUASH(IX)),IX=JX,ITOP)
 90   FORMAT(T17,A5,T27,A5,T37,A5,T47,A5,T57,A5)

      DO 10 I=1,IR
         IF (NSYM .EQ. 1) WRITE (LuOut,44) VARNAM(ISQUASH(NOPTI(I))),
     $      (X(I,J),J=JX,MIN(IC,ITOP))
         IF (NSYM .EQ. 0) WRITE (LuOut,44) VARNAM(ISQUASH(I)),
     $      (X(I,J),J=JX,MIN(IC,ITOP))
 10   CONTINUE
      JX=JX+5
      IF (JX .GT. IR) RETURN
      WRITE(LuOut,*)
      GOTO 1
 44   FORMAT(T3,A5,T13,F9.6,T23,F9.6,T33,F9.6,T43,
     &       F9.6,T53,F9.6)

      RETURN
      END
