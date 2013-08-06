      SUBROUTINE PROCESTEP(SCRATCH, HES, GRAD, SCALE, STPMAG,
     &                     STEP, STPTOL, STPMAX, IQFIX, ISTCRT,
     &                     BEGN_TRUST_RAD, EPS, QSTLST_CLIMB, TS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL TRUST_RADIUS, DELTXT_DELTX, ABSLT_LARGST, QSTLST_CLIMB
      LOGICAL TS
      DOUBLE PRECISION ONE
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)






































































































































































































































































































































































































































































































































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






c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



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
      IF (iFlags2(5) .EQ. 1) 
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
      TRUST_RADIUS = (IFLAGS(52) .EQ. 0)
      DELTXT_DELTX = (IFLAGS(52) .EQ. 1)
      ABSLT_LARGST = (IFLAGS(52) .EQ. 2)
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
            Write(6,*)
            Write(6,"(x,a,F10.5)") "The initial trust rad:", 
     &                               TRUST_RAD
            Write(6,"(x,a,a,3F10.5)") "The predicted and actual", 
     &                             " energy change: ", 
     &                               PRD_ENRG_CHNG,ACT_ENRG_CHNG
            Write(6,"(x,a,F10.5)") "Ratio of change:",
     &                              PRD_ENRG_CHNG/ACT_ENRG_CHNG
            TAU = DSQRT(DDOT(NXM6, STEP, 1, STEP, 1))
            CALL TRUST_UPDATE(SCRATCH, GRAD, HES, TRUST_RAD,
     &                        PRD_ENRG_CHNG, DXDX, EPS, NX, NXM6,
     &                        NOPT, NCYCLE, TS)
            WRITE(6,"(x,a,2F10.5)") "The new trust rad and TAU:", 
     &                             TRUST_RAD, TAU
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
