         SUBROUTINE  ERD__RYS_X_ROOTS_WEIGHTS_YUKAWA
     +
     +                    ( NT,NTGQP,
     +                      NGQP,NMOM,
     +                      TVAL,UVAL,RYSZERO,
     +                      A,B,
     +                      MOM,
     +                      DIA,OFF,
     +                      ROW1,ROW2,
     +
     +                                RTS,
     +                                WTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__RYS_X_ROOTS_WEIGHTS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : Master routine to evaluate roots and weights in the
C                interval [0,1] over the Rys weight function:
C
C
C                                       exp(-T*x)+U*(1-1/x)
C                             W   (x) = ----------------
C                              Rys           2*sqrt(x)
C
C
C                using the general Gaussian Quadrature technique
C                consisting of the following basic steps:
C
C
C                  Input:
C
C                    NT           =  # of T-values
C                    NTGQP        =  # of roots times # of T-values
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    TVAL         =  the T-values
c                    UVAL         =  the U-values
C                    RYSZERO      =  the zeroth Rys moments for all
C                                    T-values
C                    A,B          =  will contain the recurrence
C                                    coefficients for the auxillary
C                                    polynomials
C                    MOM          =  will contain the normed auxillary
C                                    polynomial modified moments
C                    DIA,OFF      =  will contain the diagonal and
C                                    offdiagonal elements of the
C                                    tridiagonal symmetric terminal
C                                    matrix
C                    ROW1,ROW2    =  will be used to evaluate the
C                                    tridiagonal elements of the
C                                    symmetric terminal matrix in an
C                                    efficient way using Sack and
C                                    Donovan's method
C
C                  Output:
C
C                    RTS          =  the roots array
C                    WTS          =  the weights array
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C                      >>>  COPYRIGHT NOTICE  <<<
C                             Program ERD
C
C    This routine is part of the integral program package ERD and as
C    such is copyrighted by N. Flocke, 2004; all rights reserved.
C    The full copyright notice is printed in the main routine of the
C    ERD package and all terms and conditions given therein are
C    incorporated here.
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INCLUDE    'erd__jacobi.inc'

         INTEGER    I,J,M,N
         INTEGER    IMAX,JMAX
         INTEGER    IP1,JP1
         INTEGER    ITER
         INTEGER    MOMMAX
         INTEGER    MXITER
         INTEGER    NGQP
         INTEGER    NMOM
         INTEGER    NRTS
         INTEGER    NT,NTGQP

         DOUBLE PRECISION   C,D,F,G,P,R,S
         DOUBLE PRECISION   BINC,SINC,STGEXP
         DOUBLE PRECISION   LIM1,LIM2,LIM3
         DOUBLE PRECISION   MOMI,MOMIM1,MOMIP1,MOMZERO
         DOUBLE PRECISION   R1
         DOUBLE PRECISION   ROOT
         DOUBLE PRECISION   SCALE
         DOUBLE PRECISION   SIGMA,THETA
         DOUBLE PRECISION   T,TINV,TINV2,TINVHF,TINVSQ,TEXP,TPOWER
         DOUBLE PRECISION   TEST1,TEST2
         DOUBLE PRECISION   ZMOM,ZINV
         DOUBLE PRECISION   ZERO,EXTREM,VSMALL,SMALL,TENTH,
     +                      HALF,ONE,ONEP5,TWO,THREE,TEN,TLIMIT

         DOUBLE PRECISION   A       (1:NMOM)
         DOUBLE PRECISION   B       (2:NMOM)
         DOUBLE PRECISION   MOM     (1:NMOM)
         DOUBLE PRECISION   DIA     (1:NGQP)
         DOUBLE PRECISION   OFF     (1:NGQP)
         DOUBLE PRECISION   ROW1    (1:NMOM)
         DOUBLE PRECISION   ROW2    (1:NMOM-1)

         DOUBLE PRECISION   RYSZERO (1:NT)
         DOUBLE PRECISION   TVAL(1:NT), UVAL(1:NT)

         DOUBLE PRECISION   RTS     (1:NTGQP)
         DOUBLE PRECISION   WTS     (1:NTGQP)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (EXTREM  =  1.0D-300)
         PARAMETER  (VSMALL  =  1.0D-16)
         PARAMETER  (SMALL   =  1.0D-6)
         PARAMETER  (TENTH   =  0.1D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (ONEP5   =  1.5D0)
         PARAMETER  (TWO     =  2.D0)
         PARAMETER  (THREE   =  3.D0)
         PARAMETER  (TEN     =  10.D0)
         PARAMETER  (TLIMIT  =  15.D0)

         PARAMETER  (MOMMAX  =  30)
         PARAMETER  (MXITER  =  30)
C
C
C------------------------------------------------------------------------|
C   Invocation of Toru Shiozaki's interpolation scheme for Rys roots     |
C       and weights using two-dimentional interpolation procedure        |
c------------------------------------------------------------------------

c         DO J=1,NT
c         UVAL(J)=1.d-7
c         write(6,*) TVAL(J),'===',UVAL(J), J
c         ENDDO

c         call root2(TVAL, UVAL, RTS, WTS, NT)

c           DO I=1,NTGQP          
c           write(6,*) RTS(I),'===',WTS(I), I 
c           ENDDO
C             ...ready!
C
C
         RETURN
         END
