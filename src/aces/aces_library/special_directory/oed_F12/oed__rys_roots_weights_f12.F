         SUBROUTINE  OED__RYS_ROOTS_WEIGHTS_F12
     +
     +                    ( NT,NTGQP,
     +                      NGQP,NMOM,
     +                      TVAL,UVAL,RYSZERO,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      A,B,
     +                      MOM,
     +                      DIA,OFF,
     +                      ROW1,ROW2,
     +
     +                             RTS,
     +                             WTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__RYS_ROOTS_WEIGHTS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__RYS_1_ROOTS_WEIGHTS
C                OED__RYS_2_ROOTS_WEIGHTS
C                OED__RYS_3_ROOTS_WEIGHTS
C                OED__RYS_4_ROOTS_WEIGHTS
C                OED__RYS_5_ROOTS_WEIGHTS
C                OED__RYS_X_ROOTS_WEIGHTS
C  DESCRIPTION : This routine calculates NGQP-point Gaussian quadrature
C                rules on [0,1] over the Rys weight functions:
C
C
C                                       exp(-T*x)
C                             W   (x) = ---------
C                              Rys      2*sqrt(x)
C
C
C                for a set of NT T-exponents. Special interpolation
C                routines are provided for low number of roots and
C                weigths (NGQP < 6). On exit, NT x NGQP = NTGQP roots
C                and weights have been produced.
C
C
C                  Input:
C
C                    NT           =  # of T-exponents
C                    NTGQP        =  # of roots times # of T-exponents
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    TVAL         =  the T-exponents
C                    RYSZERO      =  will hold the zeroth Rys moments
C                                    for all T-exponents
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    A,B          =  will contain the recurrence
C                                    coefficients for the auxillary
C                                    polynomials
C                    MOM          =  will contain the normed auxillary
C                                    polynomial modified moments
C                    DIA,OFF      =  will contain the diagonal and
C                                    offdiagonal elements of the
C                                    tridiagonal symmetric terminal
C                                    matrix
C                    ROW1,ROW2    =  first,second row intermediates.
C                                    Will be used to evaluate the
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
C                             Program OED
C
C    This routine is part of the integral program package OED and as
C    such is copyrighted by N. Flocke, 2004; all rights reserved.
C    The full copyright notice is printed in the main routine of the
C    OED package and all terms and conditions given therein are
C    incorporated here.
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    MGRID
         INTEGER    N,I
         INTEGER    NGQP
         INTEGER    NGRID
         INTEGER    NMOM
         INTEGER    NT,NTGQP
         INTEGER    TGRID

         DOUBLE PRECISION   DELTA,QQQ
         DOUBLE PRECISION   T,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION   ZERO,SIXTH,FIFTH,FOURTH,THIRD,HALF,ONE,PI

         DOUBLE PRECISION   A    (1:NMOM)
         DOUBLE PRECISION   B    (2:NMOM)
         DOUBLE PRECISION   MOM  (1:NMOM)
         DOUBLE PRECISION   DIA  (1:NGQP)
         DOUBLE PRECISION   OFF  (1:NGQP)
         DOUBLE PRECISION   ROW1 (1:NMOM)
         DOUBLE PRECISION   ROW2 (1:NMOM-1)

         DOUBLE PRECISION   RYSZERO (1:NT)
         DOUBLE PRECISION   TVAL    (1:NT)
         DOUBLE PRECISION   UVAL    (1:NT)

         DOUBLE PRECISION   RTS     (1:NTGQP)
         DOUBLE PRECISION   WTS     (1:NTGQP)

         DOUBLE PRECISION   FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (SIXTH   =  0.166666666666667D0)
         PARAMETER  (FIFTH   =  0.2D0)
         PARAMETER  (FOURTH  =  0.25D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (PI      =  3.141592653589793D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed according to the number of roots and weights
C                required.
C
c          write(6,*) '++++++++++++++++++++++++++++'
C
c------------Number of roots and weights equal to  1--------------------
         IF(NGQP.eq.1) then


c         CALL OED__RYS_1_2_ROOTS_WEIGHTS_TENNO
c     +
c     +                    ( NT,NTGQP,
c     +                      NGQP,NMOM,
c     +                      TVAL,UVAL,RYSZERO,
c     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
c     +                      A,B,
c     +                      MOM,
c     +                      DIA,OFF,
c     +                      ROW1,ROW2,
c     +
c     +                                RTS,
c     +                                WTS )
c     +

         call root1(TVAL, UVAL, RTS, WTS, NT)


         ENDIF

c------------Number of roots and weights equal to  2--------------------
         IF(NGQP.eq.2) then
c         CALL OED__RYS_1_2_ROOTS_WEIGHTS_TENNO
c     +
c     +                    ( NT,NTGQP,
c     +                      NGQP,NMOM,
c     +                      TVAL,UVAL,RYSZERO,
c     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
c     +                      A,B,
c     +                      MOM,
c     +                      DIA,OFF,
c     +                      ROW1,ROW2,
c     +
c     +                                RTS,
c     +                                WTS )
c     +


         call root2(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  3--------------------
         IF(NGQP.eq.3) then
         call root3(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  4--------------------
         IF(NGQP.eq.4) then
         call root4(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  5--------------------
         IF(NGQP.eq.5) then
         call root5(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  6--------------------
         IF(NGQP.eq.6) then
         call root6(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  7--------------------
         IF(NGQP.eq.7) then
         call root7(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  8--------------------
         IF(NGQP.eq.8) then
         call root8(TVAL, UVAL, RTS, WTS, NT)
         ENDIF




C
C             ...ready!
C
C
         RETURN
         END
