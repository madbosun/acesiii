C
C A fairly simple routine to identify symmetry unique values in the
C array RIN and put the symmetry unique elements to the begining
C of the RIN array followed by locating the elements in the RIN
C array that are being optimized and putting them in ROUT array
C in the proper (to ACES II) order.  Ajith Perera, 01/2005.
C
      SUBROUTINE SYMUNQONLY(RIN, ROUT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT

      DIMENSION RIN(NX), ROUT(NX)

      DO I = 1, NUNIQUE
         IP = IUNIQUE(I)
         Z  = RIN(IP)
         DO J = 1, NEQ(IP)
            Z = Z + RIN(IEQUIV(I,J))
         END DO
         FIAVE  = Z/(NEQ(IP)+1)
         RIN(I) = FIAVE
      END DO
      DO IUNQ = 1, NUNIQUE
         ROUT(IUNQ) = RIN(NOPTI(IUNQ))
      END DO
      RETURN
      END

