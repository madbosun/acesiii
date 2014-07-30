         SUBROUTINE  OED__NAI_DERV_INT1D_TO_A0_F12
     +
     +                    ( SHELLA,
     +                      NEXP,NGQPCEN,NGEXCEN,
     +                      NXYZA,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_INT1D_TO_A0
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative nuclear attraction integrals [A|0] or [0|A],
C                adding up the contributions from all the 1D A0 or 0A
C                integrals for all nuclear centers.
C
C                Simplified version of the general AB routine to reduce
C                loop overheads for those cases where there is at least
C                one s-shell. For comments and details see the general
C                AB routine.
C
C
C                  Input:
C
C                    SHELLA      =  shell types for contraction shell A
C                    NEXP        =  current # of exponent pairs
C                    NGQPCEN     =  product of # of gaussian quadrature
C                                   points (roots) times # of nuclear
C                                   attraction centers
C                    NGEXCEN     =  product of # of gaussian quadrature
C                                   points times # of exponent pairs
C                                   times # of nuclear attraction
C                                   centers
C                    NXYZA       =  # of cartesian monomials for shell A
C                    INT1Dx      =  all current 1D A0/0A derivative
C                                   nuclear attraction integrals for
C                                   each cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   1D A0/0A derivative integral
C                                   products
C                    SCALE       =  the NGEXCEN scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [A|0] or [0|A] derivative nuclear
C                                   attraction integrals corresponding
C                                   to all current exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C                      >>>  COPYRIGHT NOTICE  <<<
C                             Program OED
C
C    This routine is part of the integral program package OED and as
C    such is copyrighted by N. Flocke, 2004; all rights reserved.
C    The full copyright notice is printed in the main routines of the
C    OED package and all terms and conditions given therein are
C    incorporated here.
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     DIFFY,DIFFZ

         INTEGER     I,K,M,N
         INTEGER     NEXP,NGQPCEN,NGEXCEN
         INTEGER     NXYZA
         INTEGER     SHELLA
         INTEGER     XA,YA,ZA
         INTEGER     YAMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO,STGEXP,UVAL(1:NEXP),RTS(1:NGEXCEN)

         DOUBLE PRECISION  SCALE (1:NGEXCEN)
         DOUBLE PRECISION  TEMP1 (1:NGEXCEN)
         DOUBLE PRECISION  TEMP2 (1:NGEXCEN)

         DOUBLE PRECISION  BATCH (1:NEXP,1:NXYZA)

         DOUBLE PRECISION  INT1DX (1:NGEXCEN,0:SHELLA)
         DOUBLE PRECISION  INT1DY (1:NGEXCEN,0:SHELLA)
         DOUBLE PRECISION  INT1DZ (1:NGEXCEN,0:SHELLA)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...outer loop over x-contribution. No skipping of
C                x-contribution of 0-type can be done here,
C                since the 1DX integrals carry the Rys weight!
C
C
         I = 0
         DO 100 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO N = 1,NGEXCEN
               TEMP1 (N) = SCALE (N) * INT1DX (N,XA)
            END DO
C
C
C             ...inner loop over y-contribution. Skip the multiplication
C                of y-contributions, if no y-coordinate derivative
C                was formed and we have a 0-contribution, as then the
C                1DY integrals are equal to 1.
C
C
            DO 110 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGEXCEN
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGEXCEN
                      TEMP2 (N) = TEMP1 (N) * INT1DY (N,YA)
                   END DO
               END IF
C
C
C             ...skip multiplication of z-contributions, if we
C                have a 0-type and no derivations were performed
C                on the z-coordinate, as then the 1DZ integrals
C                are equal to 1. All info concerning all three
C                x,y and z contributions have been collected for all
C                exponent pairs, nuclear centers and quadrature points
C                at once. Sum up the 1D X,Y,Z integral products over
C                the nuclear centers and the quadrature points to
C                the appropriate place of the [A|0] or [0|A] batch.
C
C
               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   K = 0
                   DO M = 1,NEXP
                      SUM = ZERO
                      DO N = 1,NGQPCEN
                         SUM = SUM + TEMP2 (K+N)*(1.d0-RTS(K+N) )
                      END DO
                      K = K + NGQPCEN
                      BATCH (M,I) = SUM*(2.d0*UVAL(M)/STGEXP)
                   END DO
               ELSE
                   K = 0
                   DO M = 1,NEXP
                      SUM = ZERO
                      DO N = 1,NGQPCEN
             SUM = SUM + TEMP2 (K+N) * INT1DZ (K+N,ZA)*(1.d0-RTS(K+N) )
                      END DO
                      K = K + NGQPCEN
                      BATCH (M,I) = SUM*(2.d0*UVAL(M)/STGEXP)
                   END DO
               END IF
C
C
C             ...next z- and y-contributions.
C
C
  110       CONTINUE
C
C
C             ...next x-contribution.
C
C
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
