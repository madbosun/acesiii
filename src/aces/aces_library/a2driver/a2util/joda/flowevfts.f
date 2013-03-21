      SUBROUTINE FLOWEVFTS(SCRATCH, GRDHES, HESMOD, DIAGHES, LMBDAN, 
     &                     LMBDAP, STPMAX, MORSE, IMODE, NX, LUOUT, 
     &                     IBREAK, NOPT, QSTLST_CLIMB)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C















































































































































































































































































































































































































































































































































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
      DOUBLE PRECISION LMBDAN, LMBDAP
      LOGICAL MORSE, QSTLST_CLIMB
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), HESMOD(NOPT, NOPT),
     &          DIAGHES(NOPT, NOPT)
C
      IBREAK = 0
C
      Write(6,*)
      WRITE(6, "(a,(4F10.5))") "The step size at entry",
     &          (SCRATCH(J+NOPT), J=1, NOPT)
      WRITE(6, "(a,(4F10.5))") "The gradhes at entry",
     &          (GRDHES(J), J=1, NOPT)
      Write(6,*)

      IF (.NOT.QSTLST_CLIMB) THEN
C
      DO 10 I = 1, NOPT
C
         IF (IBREAK .EQ. 1) GOTO 10
C
         IF (I .EQ. IMODE) THEN
C
C Allow here for the possibility that the user may want break symmetry
C following a non-totally symmetry mode out of a local minimum. In this
C case Lambda(P) and eigenvalues, gradient along the cooresponding 
C eigen vectors become zero. If this occurs, the code below will force
C the geometry to follow this mode.
C
            DENOM = HESMOD(I,I) - LMBDAP


C
         ELSE
C
            DO 20 J = 1, NOPT

                IF (DABS(HESMOD(I,I)-LMBDAN) .LT. 1.0D-06) Then
                    SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
                ELSE 
                    SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-GRDHES(I)
     &                               *DIAGHES(J,I)/(HESMOD(I,I) 
     &                                - LMBDAN) 
                ENDIF
C
 20         CONTINUE
C
         ENDIF
C
 10   CONTINUE
C
C     ENDIF (.NOT.QSTLST_CLIMB)
      ENDIF
C
      Write(6,*)
      WRITE(6, "(a,(4F10.5))") "The step size before Morse",
     &          (SCRATCH(J+NOPT), J=1, NOPT)
      Write(6,*)

      IF (MORSE) THEN
         IF (iFlags2(5) .ge. 3) THEN
            CALL DOMORSEXYZ_RIC(SCRATCH, NOPT, NX, LUOUT)
         ELSE IF (iFlags2(5) .eq. 1) THEN
           CALL DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
         ENDIF 
      ENDIF
C
C Add in part of step which goes along the reaction coordinate.
C
      DO 30 J = 1, NOPT
         IF (QSTLST_CLIMB) DENOM = HESMOD(IMODE, IMODE) - LMBDAP
         
         IF (DABS(DENOM) .LT. 1.0D-06) THEN
            SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
         ELSE
            SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-GRDHES(IMODE)
     &                        *DIAGHES(J,IMODE)/DENOM
         ENDIF
C     
 30   CONTINUE
C
      Write(6,*)
      WRITE(6,"(a,(4F10.5))") "The unscaled step size", 
     &            (SCRATCH(J+NOPT), J=1, NOPT)
      Write(6,*)
C
      RETURN
      END
