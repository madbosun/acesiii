      SUBROUTINE NUCREP(ZNREP,IERR)
C
C DETERMINES THE INTERNUCLEAR REPULSION ENERGY.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Maximum number of atoms currently allowed
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


C
      IQ(I)=3*I-2
      Z=0.D0
      DO 10 I=1,NATOMS
      DO 10 J=I+1,NATOMS
       BOT=DIST(Q(IQ(I)),Q(IQ(J)))
       IF(BOT.LT.1.D-5)THEN
        IERR=1
        GOTO 10
       ELSE
        IF(IATNUM(I).NE.110.AND.IATNUM(J).NE.110)THEN
c        IGNORE GHOST ATOMS
         TOP=FLOAT(IATNUM(I)*IATNUM(J))
         Z=Z+(TOP/BOT)
        ENDIF
      ENDIF
10    CONTINUE
      ZNREP=Z
      RETURN
      END
