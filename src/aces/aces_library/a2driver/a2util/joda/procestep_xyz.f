      SUBROUTINE PROCESTEP_XYZ(SCRATCH, AMATRX, TOTREDNCO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL CONVERGED, CHANGED, NO_ITER
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG
      DOUBLE PRECISION ORIENT(3,3)
      PARAMETER (MAXITER = 100)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)


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















































































































































































































































































































































































































































































































c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
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
      PARAMETER (EPSILON = 1.0D-5)
      CHARACTER*80 HOLDQ
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION SCRATCH(9*NATOMS*NATOMS), IQFIX(3*MXATMS, 3*MXATMS),
     &          AMATRX(3*NATOMS*TOTREDNCO), STATS(6), 
     &          REDUNCO(MAXREDUNCO), SCRATCH_LCL(4*MAXREDUNCO+6*
     &          MXATMS+15*MXATMS),
     &          TMP1(MAXREDUNCO*MAXREDUNCO),TMP2(MAXREDUNCO)
C
      DATA BTOA /0.529177249D0/

      DO I = 1, NUNIQUE
         DO J = 1, NEQ(IUNIQUE(I))
            IQFIX(IUNIQUE(I), J) = IEQUIV(I,J)
         ENDDO
      ENDDO
C
      DO J = 1, NOPT
C
         SCRATCH_LCL(NOPTI(J)) = SCRATCH(J + NOPT)
         SCRATCH_LCL(NOPT + NOPTI(J)) = SCRATCH(J + NOPT)
C
         DO K = 1, NEQ(NOPTI(J))
            SCRATCH_LCL(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
           SCRATCH_LCL(NOPT + IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)  
         ENDDO
      ENDDO
C
      Write(6,*)
      WRITE(6,*) "The optimization parameters:nx, nopt, totrednco nxm6",
     &            nx, totrednco, nopt, nxm6
      Write(6,*)
      WRITE(6,*) "The delta q0"
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(I), I = 1, TOTREDNCO)
      Write(6,*) "The The rnew (q1)"
      WRITE(6,10) (R(I), I = 1, TOTREDNCO)
   10 Format (5(1X,F10.6))
      Write(6,*)
C
      IOFF4q0  = 3*NATOMS + 3*TOTREDNCO
      IOFF4X0  = IOFF4q0  + TOTREDNCO
C
      ISCRATCH = 4*TOTREDNCO + 6*NATOMS 
      ITMPQ    = ISCRATCH    + 9*NATOMS
      INEWQ    = ITMPQ       + 3*NATOMS
      ILAST    = INEWQ       + 3*NATOMS
C
      CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL, 1, R, 1)
      CALL DCOPY(TOTREDNCO, R, 1, SCRATCH_LCL(2*TOTREDNCO + 1), 1)
      CALL DCOPY(TOTREDNCO, SCRATCH_LCL, 1, SCRATCH_LCL(IOFF4q0 + 1),
     &           1)
      CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(IOFF4X0 + 1), 1)

C#ifdef _TEST_RED
C This requires defining TMP1 and TMP2 arrays.
      CALL GETREC(20,'JOBARC','NUMREDCO', 1, NULLEVAL)
      CALL GETREC(20,'JOBARC','REDEVECS',NULLEVAL*TOTREDNCO*IINTFP,
     &            TMP1)
      CALL XGEMM('T', 'N', NULLEVAL, 1, TOTREDNCO, 1.0D0,
     &               TMP1, TOTREDNCO, SCRATCH_LCL, TOTREDNCO,
     &               0.0D0, TMP2, NULLEVAL)
      Write(6,*) "Check for the redundencies in deltaq0"
      Write(6,"(6F13.7)") (tmp2(i), i=1, NULLEVAL)
C#endif
C
C     
      CONVERGED = .FALSE. 
        CHANGED = .FALSE.
        NO_ITER = .FALSE.
          NITER = 0 
C
      DO WHILE (.NOT. CONVERGED .AND. NITER .LE. MAXITER .AND.
     &          .NOT. CHANGED .AND. .NOT. NO_ITER) 
C
         CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(3*TOTREDNCO + 1), 1)
         Write(6,*)
CSSS         Write(6,*) "The A matrix"
CSSS         CALL OUTPUT(AMATRX, 1, 3*NATOMS, 1, TOTREDNCO, 3*NATOMS,
CSSS     &               TOTREDNCO, 1) 
         WRITE(6,*) "The starting Cartesian during iteration; x_K"
         CALL OUTPUT(Q, 1, 3, 1, Natoms,3, Natoms, 1)
         Write(6,*)
         Write(6,*) "The starting internal increments:delta q_k"
         Write(6,*) 
         WRITE(6,10) (SCRATCH_LCL(I), I=1, TOTREDNCO)
         Write(6,*)
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &               AMATRX, 3*NATOMS, SCRATCH_LCL, TOTREDNCO,
     &               1.0D0, Q, 3*NATOMS)
C
CSSS         CALL GETREC(20,'JOBARC','ORIENT2 ',9*IINTFP,ORIENT)
CSSS         CALL XGEMM('T', 'N', 3, NATOMS, 3, 1.0D0, ORIENT, 3, Q,
CSSS     &               3, 0.0D0, SCRATCH_LCL(INEWQ), 3)
CSSS         CALL  DCOPY(3*NATOMS, SCRATCH_LCL(INEWQ), 1, Q, 1)
C
         Print*, "The Cartesian after the update: x_k+1"
         CALL OUTPUT(Q, 1, 3, 1, Natoms,3, Natoms, 1)
         CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNOFBND,
     &                    NW_TOTREDNCO, TOTNOFANG, .FALSE.)
C
         IF (IFLAGS(60).EQ. 2) THEN
            CALL GETREC(20,'JOBARC','PLSMINSP',NW_TOTREDNCO,NCON)  
            DO ICOORD = 1, NW_TOTREDNCO
               IF (NCON(ICOORD).EQ.1) REDUNCO(ICOORD) = -REDUNCO(ICOORD)
            END DO
         END IF
C
C In rare cases, it is possible that the # of RIC coordinates
C can change. If that happens leave the iterative procedure and
C and use the linear update.
C
         IF (NW_TOTREDNCO .GT. TOTREDNCO) CHANGED = .TRUE.
         CALL CHECK4_180(SCRATCH_LCL(2*TOTREDNCO + 1), SCRATCH_LCL,
     &                   TOTNOFBND, TOTNOFANG, TOTREDNCO, NO_ITER)
C        
C
         CALL DCOPY(TOTREDNCO, REDUNCO, 1, R, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL(2*TOTREDNCO + 1),
     &              1, R, 1)
         CALL DCOPY(TOTREDNCO, SCRATCH_LCL(IOFF4q0 + 1), 1, 
     &              SCRATCH_LCL, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, R, 1, SCRATCH_LCL, 1)
C
      WRITE(6,*) "The delta delta q_k" 
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(I), I=1, TOTREDNCO)
      Write(6,*)
         CALL REMOVE_360(SCRATCH_LCL, TOTNOFBND, TOTNOFANG, 
     &                   TOTREDNCO)
C
         DO IREDUNCO = 1, TOTREDNCO
            IF (DABS(SCRATCH_LCL(IREDUNCO)) .LT. 1.0D-10) SCRATCH_LCL
     &          (IREDUNCO) = 0.0D0
         ENDDo
   
         IF (NITER .GE. 0) CALL VSTAT(SCRATCH_LCL, STATS, TOTREDNCO)
         IF (NITER .EQ. 0) THEN
             RMSOF_DDQ1 = STATS(5)
         ELSE
C
CSSS             IF (STATS(5) .GE. RMSOF_DDQ1) NO_ITER = .TRUE.
             IF (STATS(5) .GE. RMSOF_DDQ1) NO_ITER = .FALSE.
         ENDIF
C
C#ifdef _TEST_RED
C This requires defining TMP1 and TMP2 arrays.
      CALL XGEMM('T', 'N', NULLEVAL, 1, TOTREDNCO, 1.0D0,
     &               TMP1, TOTREDNCO, SCRATCH_LCL, TOTREDNCO,
     &               0.0D0, TMP2, NULLEVAL)
      Write(6,*)
      Write(6,*) "Check for the redundencies in delta delta q_K"
      Write(6,"(6F13.7)") (tmp2(i), i=1, NULLEVAL)
C#endif
C
         CALL DAXPY(3*NATOMS, -1.0D0, Q, 1, SCRATCH_LCL
     &             (3*TOTREDNCO + 1), 1)
         CALL VSTAT(SCRATCH_LCL(3*TOTREDNCO+1), STATS, 3*NATOMS)
C
         IF (STATS(5)  .LE. EPSILON) CONVERGED = .TRUE.
         NITER = NITER + 1 

      ENDDO 
      IF (CONVERGED) Write(6, "(a,a,I2,a)") "Iterative update ",
     &                         "converged in ", NITER, " iterations."
C
C If the iterative procedure fail to converge in maximum iterations
C alowed, then do a simple linear update.
C
C#endif
C
      IF (.NOT. CONVERGED) THEN
         Write(6, "(a,a)") "Iterative update did not converged.",
     &                     " linear update is used. "
         CALL DCOPY(3*NATOMS, SCRATCH_LCL(IOFF4X0 +1), 1, Q, 1)
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &              AMATRX, 3*NATOMS, SCRATCH_LCL(IOFF4q0 + 1), 
     &              TOTREDNCO, 1.0D0, Q, 3*NATOMS)
      ENDIF
C
C Lets screen the Cartesian coordiantes and filter out changes
C beyond 7th decimal.
C 
      DO IATOM = 1, 3*NATOMS
         WRITE(HOLDQ, "(F12.7)") Q(IATOM)
         READ(HOLDQ,  "(F12.7)") Q(IATOM)
      ENDDO 

      Write(6,*)
      WRITE(6,*) "The new Cartesian Coord., x+delta x0"
      WRITE(6,"(3F13.9)") (Q(I), I = 1, 3*NATOMS)
C
      CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNOFBND,
     &                 NW_TOTREDNCO, TOTNOFANG, .TRUE.)

C
C Store the Cartesian coordinates with a different record label
C from COORD since finite difference calculations use COORD
C for the grid points, not the current point.
C

      CALL PUTREC(1,'JOBARC','COORD_OP',IINTFP*3*NATOMS,Q)
C
      RETURN
      END
