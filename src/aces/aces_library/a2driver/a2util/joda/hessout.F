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
      
CSS     IF(NSYM.EQ.1)
CSS     $   WRITE(LuOut,90)(VARNAM(ISQUASH(NOPTI(IX))),IX=JX,ITOP)
CSS      IF(NSYM.EQ.0)WRITE(LuOut,90)(VARNAM(ISQUASH(IX)),IX=JX,ITOP)
CSS 90   FORMAT(T15,A5,T25,A5,T35,A5,T45,A5,T55,A5)
CSS
      IF (iFlags2(h_IFLAGS2_geom_opt) .eq. 2) THEN
         CALL OUTPUT(X, 1, NX, 1, NX, NX, NX, 1)
      ELSE
         CALL OUTPUT(X, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
      ENDIF
       
#ifdef _NOSKIP
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
 44   FORMAT(T3,A5,T12,F12.6,T22,F12.6,T32,F12.6,T42,
     &       F12.6,T52,F12.6)
#endif 
      RETURN
      END
