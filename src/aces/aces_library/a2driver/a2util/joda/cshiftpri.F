C
      SUBROUTINE CSHIFTPRI(CSHIFT, INDEX)
C
C Print out the total shielding tensor and orientation of the shielding
C tensor. Coded by Ajith 04/96.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION CSHIFT (3, 3)
C
      CHARACTER*5 COORD(3)
C
#include "mxatms.par"
#include "cbchar.com"
#include "coord.com"
C
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
C
      DATA COORD/'x    ','y    ','z    '/
C
      WRITE (6, 10)
C
      DO 5 I = 1, 3
         WRITE(6, 20) ZSYM(INDEX), COORD(I), (CSHIFT(I,J),J=1,3)
    5 CONTINUE
C
 10   FORMAT(/)
 20   FORMAT(15X, A2, 1X, A2, 3F12.6)
C
      RETURN
      END
