      SUBROUTINE NEWRAPH(SCRATCH, GRDHES, HESMOD, DIAGHES, MORSE, 
     &                   NOPT, NX, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
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






























































































































































































































































































































































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C
      LOGICAL MORSE, XYZIN, NWFINDIF
C
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), DIAGHES(NOPT, NOPT),
     &          HESMOD(NOPT, NOPT)
C     
C Do the Normal Newton-Raphson update
C
      Write(6,*) 
      Write(6, "(1x,a,I5)") "The number of degs. of freed. 
     &          at NR:", NOPT
      Write(6,*)
      DO 10 I = 1, NOPT
         DO 20 J = 1, NOPT
            IF (HESMOD(I,I) .EQ. 0.0D0) Then
                    SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
            ELSE
C
                SCRATCH(J+NOPT)=SCRATCH(J+NOPT)-GRDHES(I)*
     &                          DIAGHES(J,I)/HESMOD(I,I) 
            ENDIF
 20      CONTINUE
 10   CONTINUE
C
      Write(6,*) "The unmolested step"
      Write(6, "(8F10.4)") (SCRATCH(J+NOPT), J=1, NOPT)
C
C Notice that there is no Morse adjustments can be made for
C pure Cartesian OLPtimization (no connections), Also, no 
C Morse scaling can be done for redundent internals.
C Ajith Perera,08/2008, 2011.
C
      IF (MORSE) THEN
         IF (iFlags2(5) .ge. 3) THEN
CSSS           CALL DOMORSEXYZ_RIC(SCRATCH, NOPT, NX, LUOUT)
         ELSE IF (iFlags2(5) .eq. 1) THEN
           CALL DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
         ENDIF  
      ENDIF
C
      RETURN
      END
