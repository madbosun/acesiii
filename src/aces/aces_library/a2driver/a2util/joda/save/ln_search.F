      SUBROUTINE LINSRC(ICRD,MCRD,GRAD,GRADOL,STEP,STEPOL,
     &     TMPVEC,TMPVC2,ACTIVE,EMOD)
C
C     This routine calculates the step that should be taken to obtain
C     the next geometry.
C
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0)
      LOGICAL ACTIVE
      DIMENSION GRAD(MCRD), GRADOL(MCRD), STEP(MCRD), STEPOL(MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      ACTIVE = .TRUE.
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(STEP,MCRD)
      STNM = SQRT(DDOT(ICRD,STEPOL,1,STEPOL,1))
      DO 10 I = 1, ICRD
         TMPVEC(I) = -STEPOL(I)/STNM
 10   CONTINUE
      GRAD0 = DDOT(ICRD,GRADOL,1,TMPVEC,1)
      GRAD1 = DDOT(ICRD,GRAD,1,TMPVEC,1)
      IF ((ABS(GRAD0) .LT. 1.0D-6) .OR. (ABS(GRAD1) .LT. 1.0D-6)) THEN
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Too close to minimum, line search skipped.'
         END IF
         ACTIVE = .FALSE.
         RETURN
      END IF
C
C     The energy and gradient from the last and current point, are
C     fitted to a quartic polynomial.
C
      CA = ABS(5.0D0*GRAD0 - 2.0D0*GRAD1 + 3.0D0*ERGOLD - 3.0D0*ENERGY)
      CB = GRAD1 + GRAD0 + 2.0D0*ERGOLD - 2.0D0*ENERGY - 2.0D0*CA
      CC = ENERGY - ERGOLD - GRAD0 - CA - CB
      CD = GRAD0
      CE = ERGOLD
C
C     The line search methdod is only used if the minimum lies
C     between the two points.
C
      IF ((GRAD0 .LT. D0) .AND. (GRAD1 .GT. D0)) THEN
         THRG = 1.0D-5*MIN(ABS(GRAD0),ABS(GRAD1))
         CRDA = D0
         CRDB = 1.0D0
         GRDA = GRAD0
         GRDB = GRAD1
         ISAFE = 0
 15      CONTINUE
         ISAFE = ISAFE + 1
         IF (ISAFE .GE. 200) THEN
            ACTIVE = .FALSE.
            IF (IPRINT .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Line search failed, ignoring.'
            END IF
            RETURN
         END IF
         CRDC = CRDA + (CRDB-CRDA)*
     &        MAX(0.1D0,MIN(0.9D0,ABS(GRDA/GRDB)*0.5D0))
         GRDC = ((4.0D0*CA*CRDC + 3.0D0*CB)*CRDC + 2.0D0*CC)*CRDC + CD
         IF (ABS(GRDC) .GT. THRG) THEN
            IF (GRDC .GT. D0) THEN
               CRDB = CRDC
               GRDB = GRDC
            ELSE
               CRDA = CRDC
               GRDA = GRDC
            END IF
            GOTO 15
         END IF
C
C     If the line search ends up almost back at the previous point, we do
C     not trust it and simply discard it.
C
         IF (CRDC .LT. 0.15D0) CRDC = 1.D0
C
         DO 20 I = 1, ICRD
            STEP(I) = STEPOL(I)*(1.0D0-CRDC)
            GRAD(I) = GRADOL(I) + (GRAD(I)-GRADOL(I))*CRDC
 20      CONTINUE
         EMOD = (((CA*CRDC + CB)*CRDC + CC)*CRDC + CD)*CRDC + CE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Interpolated step',-1)
            CALL OUTPUT(STEP,1,1,1,ICRD,1,MCRD,1,LUPRI)
            CALL HEADER('Interpolated gradient',-1)
            CALL OUTPUT(GRAD,1,1,1,ICRD,1,MCRD,1,LUPRI)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Interpolated energy: ',EMOD
         END IF
      ELSE
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Line search skipped.'
         END IF
         ACTIVE = .FALSE.
         RETURN
      END IF
      RETURN
      END

