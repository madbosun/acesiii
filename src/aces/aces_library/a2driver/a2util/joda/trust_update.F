
      SUBROUTINE TRUST_UPDATE(SCRATCH, GRADINT, HESIAN, TRUST_RAD,
     &                        PRD_ENRG_CHNG, DXDX, EPS, NX, NXM6, NOPT,
     &                        NCYCLE, TS)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      DOUBLE PRECISION MIN_RMIN,  MIN_RMAX, MIN_ENRG_CHNG
      LOGICAL TS

#include "machsp.com"

      DIMENSION SCRATCH(NX*NX), GRADINT(NOPT), HESIAN(NOPT, NOPT)

      DATA HALF /0.50D0/, TS_RMIN /0.70D0/, TS_RMAX /0.85D0/,
     &     MIN_RMIN /0.25D0/, MIN_RMAX /0.75d0/, TWO /2.0D0/,
     &     SCALE_FCTR /2.0D0/, ZERO /0.0D0/, ONE /1.0D0/,
     &     MIN_ENRG_CHNG /1.0D-6/, STP_MIN /0.05D0/,
     &     STP_MAX /0.60D0/, THREE /3.0D0/

      CALL GETREC(1,'JOBARC','TOTENERG',IINTFP,CURR_ENRG)
      CALL GETREC(1,'JOBARC','OLDENERG',IINTFP,PREV_ENRG)

      IF (NCYCLE .GE. 2) THEN
         CALL GETREC(1,'JOBARC','T_RADIUS',IINTFP,TAU)
         TRUST_RAD = TAU
      END IF

CSSS      WRITE(6,*) "The energy change", CURR_ENRG, PREV_ENRG,
CSSS     &            PRD_ENRG_CHNG, STP_MAX , TRUST_RAD
C
C The actual energy change
C
      ACT_ENRG_CHNG = CURR_ENRG - PREV_ENRG
C
CSSS      WRITE(6,*) "The current energy change", ACT_ENRG_CHNG,
CSSS     &            PRD_ENRG_CHNG
C
      IF (DABS(PRD_ENRG_CHNG).GT.EPS) THEN
         CURR_TRUST_RAD = ACT_ENRG_CHNG/PRD_ENRG_CHNG
         IF (TS) THEN
            IF (CURR_TRUST_RAD .LE. TS_RMIN .OR. CURR_TRUST_RAD .GE.
     &         (TWO - TS_RMIN)) TRUST_RAD = TRUST_RAD/SCALE_FCTR
            IF (CURR_TRUST_RAD .GE. TS_RMAX .OR. CURR_TRUST_RAD .LE.
     &         (TWO - TS_RMAX)) TRUST_RAD = TRUST_RAD*
     &                                      DSQRT(SCALE_FCTR)
         ELSE
CSSS            Write(6,*) "Updating trust radius", SCALE_FCTR, MIN_RMIN,
CSSS     &        CURR_TRUST_RAD, TRUST_RAD
            IF (CURR_TRUST_RAD .LE. MIN_RMIN .OR.
     &          CURR_TRUST_RAD .GT. THREE) THEN
               TRUST_RAD = TRUST_RAD/SCALE_FCTR
            ELSE IF (CURR_TRUST_RAD.GE.MIN_RMAX.AND.
     &               (DABS(ACT_ENRG_CHNG).GE.MIN_ENRG_CHNG.OR.
     &                DABS(PRD_ENRG_CHNG).GE.MIN_ENRG_CHNG)) THEN
               TMP_TRUST_RAD = TRUST_RAD*DSQRT(SCALE_FCTR)
               TMP_TRUST_RAD = MIN(TMP_TRUST_RAD,THREE*DXDX)
               TRUST_RAD = MAX(TRUST_RAD,TMP_TRUST_RAD)
            END IF
         END IF

C###         TRUST_RAD = MAX(TRUST_RAD, STP_MIN)
C###         TRUST_RAD = MIN(TRUST_RAD, STP_MAX)

       ELSE
C
C Make sure that unrealisticaly large or small steps are not taken
C
          TRUST_RAD = MAX(TRUST_RAD, STP_MIN)
          TRUST_RAD = MIN(TRUST_RAD, STP_MAX)
       END IF
C
C Calculate the "new" predicted energy change (Grad DeltaX + 1/2 DelatX H DeltaX)
C
CSSS      CALL XGEMM('N', 'N', NOPT, 1, NOPT, ONE, HESIAN, NOPT,
CSSS     &           SCRATCH(NOPT + 1), NOPT, ZERO, SCRATCH(2*NOPT +1),
CSSS     &           NOPT)
C
C The "new" predicted energy change.
C
CSSS      PRD_ENRG_CHNG = DDOT(NOPT, GRADINT, 1, SCRATCH(NOPT + 1), 1) +
CSSS     &                HALF*DDOT(NOPT, SCRATCH(2*NOPT + 1), 1,
CSSS     &                SCRATCH(NOPT + 1), 1)
C
CSSS       Write(6,*) "Leaving TRUST_RAD and Energy", TRUST_RAD,
CSSS     &             PRD_ENRG_CHNG
C
CSSS       TAU = DSQRT(DABS(DDOT(NXM6, SCRATCH(NOPT + 1), 1, SCRATCH(NOPT
CSSS     &             + 1), 1)))
C
C What we should save for the next cycle is something not clear to me.
C Nor It is clear in the papers. At this point what I have it consistent
C with MOPAC. Gamess uses the TAU of the current iteration.

       CALL PUTREC(1,'JOBARC','T_RADIUS',IINTFP,TRUST_RAD)

       RETURN
       END

