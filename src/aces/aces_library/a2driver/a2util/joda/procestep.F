      SUBROUTINE PROCESTEP(SCRATCH, HES, GRAD, SCALE, STPMAG,
     &                     STEP, STPTOL, STPMAX, IQFIX, ISTCRT,
     &                     BEGN_TRUST_RAD, EPS, QSTLST_CLIMB, TS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL TRUST_RADIUS, DELTXT_DELTX, ABSLT_LARGST, QSTLST_CLIMB
      LOGICAL TS
      DOUBLE PRECISION ONE
C
#include "mxatms.par"
#include "flags.h"
#include "jodaflags.com"
#include "machsp.com"
#include "coord.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
C
      DIMENSION SCRATCH(NX*NX), IQFIX(3*MXATMS, 3*MXATMS), AV(6),
     &          STEP(NXM6), GRAD(NOPT), HES(NOPT, NOPT), 
     &          ITIE(3*MXATMS)
      DATA ONE /1.0D0/, THREE /3.0D0/
C
C Convert to Internal from symmetry coordinates.
C
      DO 10 I = 1, NOPT
         SCRATCH(I+NOPT)=SCRATCH(I+NOPT)/DSQRT(DFLOAT(NEQ(NOPTI(I))+1))
 10   CONTINUE
C
C Filter step - Remove very small steps which might lead to broken
C symmetry.
C (STPTOL is set 10^(-12), note that STPTOL has nothing to do with STPMAX).
C
      CALL FILTER(SCRATCH(NOPT + 1), NOPT, STPTOL)
      IF (iFlags2(h_IFLAGS2_geom_opt) .EQ. 1) 
     &    CALL TIEDIHED(SCRATCH(NOPT+1), NOPT, ITIE)
C
C Step size control:
C    The old step size control algorithms were based on the largest
C absolute step size or the norm of the displacement vector (bearing
C some resemblance to the trust radius method). Improvements have been
C made to step size control based on the trust radius (Fletcher and Reeves,
C Comput. J. 13, 185 (1972)). The current implementation is based on
C P. Y. Ayala and B. S. Schlegel, JCP, 107, 375, 1997).
C
      TRUST_RADIUS = (IFLAGS(h_IFLAGS_stp_siz_ctl) .EQ. 0)
      DELTXT_DELTX = (IFLAGS(h_IFLAGS_stp_siz_ctl) .EQ. 1)
      ABSLT_LARGST = (IFLAGS(h_IFLAGS_stp_siz_ctl) .EQ. 2)
C
C Unpack the scratch array and regenerate the full STEP vector. Also
C unpack the step vector to avoid any problems with Hessian updates
C in the next cycle (this is repeated after the scaling is done below).
C I have a problem with employing symmetry-packed vectors for step size
C control as was done in the past. Only for trust radius updates am I
C undoing this. I will leave what was done in the past the way they were
C in order not to change too much. Ajith Perera, 11/2004.
C
      IF (TRUST_RADIUS) THEN
         DO J = 1, NOPT
            STEP(NOPTI(J)) = SCRATCH(J + NOPT)
            DO K=1, NEQ(NOPTI(J))
               STEP(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
            END DO
         END DO
      END IF
C
#ifdef _DEBUG_LVLM1
      Write(6,"(a)") "@-PROCESTEP, The step after symm. expansion"
      Write(6,"(F10.5)")  (STEP(J), J=1, NOPT)
#endif
C
      IF (NCYCLE.GE.2) THEN

         IF (TRUST_RADIUS) THEN

            CALL GETREC(1,'JOBARC','T_RADIUS',IINTFP,TRUST_RAD)
            CALL GETREC(1,'JOBARC','PRDENCHN',IINTFP,PRD_ENRG_CHNG)
            CALL GETREC(1,'JOBARC','PRVIUTAU',IINTFP,DXDX)
            CALL GETREC(20, 'JOBARC', 'TOTENERG', IINTFP, CURR_ENRG)
            CALL GETREC(20, 'JOBARC', 'OLDENERG', IINTFP, PREV_ENRG)
            ACT_ENRG_CHNG = CURR_ENRG - PREV_ENRG 
C
#ifdef _DEBUG_LVL0
            Write(6,*)
            Write(6,"(x,a,F10.5)") "The initial trust rad:", 
     &                               TRUST_RAD
            Write(6,"(x,a,a,3F10.5)") "The predicted and actual", 
     &                             " energy change: ", 
     &                               PRD_ENRG_CHNG,ACT_ENRG_CHNG
            Write(6,"(x,a,F10.5)") "Ratio of change:",
     &                              PRD_ENRG_CHNG/ACT_ENRG_CHNG
#endif
            TAU = DSQRT(DDOT(NXM6, STEP, 1, STEP, 1))
            CALL TRUST_UPDATE(SCRATCH, GRAD, HES, TRUST_RAD,
     &                        PRD_ENRG_CHNG, DXDX, EPS, NX, NXM6,
     &                        NOPT, NCYCLE, TS)
#ifdef _DEBUG_LVL0
            WRITE(6,"(x,a,2F10.5)") "The new trust rad and TAU:", 
     &                             TRUST_RAD, TAU
#endif
            IF (TAU.GT.TRUST_RAD) THEN
                SCALE = TRUST_RAD/TAU
            ELSE
                SCALE = ONE
            END IF

         ELSE IF (DELTXT_DELTX) THEN

            dtmp = xdot(NOPT,SCRATCH(NOPT+1),1,SCRATCH(NOPT+1),1)
            STPMAG = DSQRT(dtmp)
            SCALE = MIN(STPMAX/MAX(STPMAG,1.0d-7), 1.0D0)

         ELSE IF (ABSLT_LARGST) THEN

            CALL VSTAT(SCRATCH(NOPT + 1), AV(1), NOPT)
            SCALE = MIN(STPMAX/MAX(AV(1), 1.0D-7), 1.0D0)
            STPMAG=MAX(1.0D-7, AV(1))

         END IF

C     ELSE IF (NCYCLE.LT.2) THEN
      ELSE

         IF (TRUST_RADIUS) THEN

            IF (NCYCLE.EQ.1) THEN
               CALL PUTREC(1,'JOBARC','T_RADIUS',IINTFP,BEGN_TRUST_RAD)
            END IF
            TAU = DSQRT(DDOT(NXM6, STEP, 1, STEP, 1))
            IF (TAU.GT.STPMAX) THEN
               SCALE = BEGN_TRUST_RAD/TAU
            ELSE
               SCALE = ONE
            END IF
            STPMAG = TAU
         ELSE IF (DELTXT_DELTX) THEN

            dtmp = xdot(NOPT,SCRATCH(NOPT+1),1,SCRATCH(NOPT+1),1)
            STPMAG = DSQRT(dtmp)
            SCALE = MIN(STPMAX/MAX(STPMAG,1.0d-7), 1.0D0)

         ELSE IF (ABSLT_LARGST) THEN

            CALL VSTAT(SCRATCH(NOPT + 1), AV(1), NOPT)
            SCALE = MIN(STPMAX/MAX(AV(1), 1.0D-7), 1.0D0)
            STPMAG=MAX(1.0D-7, AV(1))

         END IF

      END IF
C
C Update the R vector
C
      IF (SCALE .NE. 1.0D0) CALL xscal(NOPT,SCALE,SCRATCH(NOPT+1),1)
C
      CALL PUTREC(1,'JOBARC','OLDGEOMT',NXM6*IINTFP,SCRATCH(NOPT+1))
C
      CALL VADD(SCRATCH(1), SCRATCH(1), SCRATCH(NOPT+1), NOPT, 1.0D0)
C
C Unpak the scratch array and regenerate the full R vector. Also
C unpak the step vector to avoid any problems with Hessian updates
C in the next cycle. 
C
      DO 400 J = 1, NOPT
C
         R(NOPTI(J)) = SCRATCH(J)
         STEP(NOPTI(J)) = SCRATCH(J + NOPT)
C
         DO 410 K=1, NEQ(NOPTI(J))
C
            STEP(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
            R(IQFIX(NOPTI(J), K)) = SCRATCH(J)
C
 410     CONTINUE
 400  CONTINUE
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      WRITE(6,*) "The scaling factor", scale
      Write(6,*) "The R vector to check I can have the old R"
      Write(6,"(5F13.7)"), (R(I), I =1, NXM6)
#endif
C
C
C Calculate the predicted energy change (Grad DeltaX + 1/2 DelatX H DeltaX)
C
      IF (TRUST_RADIUS) THEN

         CALL XGEMM('N','N',NOPT,1,NOPT,
     &              ONE, HES,              NOPT,
     &                   SCRATCH(1+NOPT),  NOPT,
     &              ZERO,SCRATCH(1+NOPT*2),NOPT)
C
C The predicted energy change.
C
         PRD_ENRG_CHNG = DDOT(NOPT, GRAD, 1, SCRATCH(NOPT + 1), 1) +
     &                   0.5d0*DDOT(NOPT, SCRATCH(2*NOPT + 1), 1,
     &                   SCRATCH(NOPT + 1), 1)
         TAU = DSQRT(DDOT(NXM6, STEP, 1, STEP, 1))
         STPMAG = TAU
CSSS         Write(6,*) "The abs step size in cycles :", TAU, STPMAX
         CALL PUTREC(1,'JOBARC','PRDENCHN',IINTFP,PRD_ENRG_CHNG)
         CALL PUTREC(1,'JOBARC','PRVIUTAU',IINTFP,TAU)

      END IF

      RETURN
      END
