C
C A fairly simple routine to identify symmetry unique values in the
C array RIN and put the symmetry unique elements to the begining
C of the RIN array followed by locating the elements in the RIN
C array that are being optimized and putting them in ROUT array
C in the proper (to ACES II) order.  Ajith Perera, 01/2005.
C
      SUBROUTINE SYMUNQONLY(RIN, ROUT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end



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

