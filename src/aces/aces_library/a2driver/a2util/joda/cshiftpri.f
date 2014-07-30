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
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


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
