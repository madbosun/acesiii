         SUBROUTINE  ERD__RYS_ROOTS_WEIGHTS_YUKAWA
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
C  OPERATION   : ERD__RYS_ROOTS_WEIGHTS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__RYS_1_ROOTS_WEIGHTS
C                ERD__RYS_2_ROOTS_WEIGHTS
C                ERD__RYS_3_ROOTS_WEIGHTS
C                ERD__RYS_4_ROOTS_WEIGHTS
C                ERD__RYS_5_ROOTS_WEIGHTS
C                ERD__RYS_X_ROOTS_WEIGHTS
C  DESCRIPTION : This routine calculates NGQP-point Gaussian quadrature
C                rules on [0,1] over the Rys weight functions:
C
C
C                                       exp(-T*x)+U*(1-1/x)
C                             W   (x) = -------------------
C                              Rys          2*sqrt(x)*x
C
C
c                          where U = (STGEXP^2)/4*rho
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
c                    UVAL         =  the U-exponents 
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

         INTEGER    MGRID
         INTEGER    N,I,K,L
         INTEGER    NGQP
         INTEGER    NGRID
         INTEGER    NMOM
         INTEGER    NT,NTGQP
         INTEGER    TGRID

         DOUBLE PRECISION   DELTA,STGEXP
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
c         write(6,*) 'Number of roots and weights = ',NGQP
c         write(6,*) 'Total number of roots and weights for all T and U',
c     & NTGQP 
C
c         DO I=1,NT
c         UVAL(I) = 1.d-7
c         ENDDO

c------------Number of roots and weights equal to  1--------------------
         IF(NGQP.eq.1) then

c        CALL ERD__RYS_1_2_ROOTS_WEIGHTS_TENNO
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
c         CALL ERD__RYS_1_2_ROOTS_WEIGHTS_TENNO
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
c
c          write(6,*) 'Ten-no roots and weights'



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

c------------Number of roots and weights equal to  9--------------------
         IF(NGQP.eq.9) then
         call root9(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  10--------------------
         IF(NGQP.eq.10) then
         call root10(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  11--------------------
         IF(NGQP.eq.11) then
         call root11(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  12--------------------
         IF(NGQP.eq.12) then
         call root12(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  13--------------------
         IF(NGQP.eq.13) then
         call root13(TVAL, UVAL, RTS, WTS, NT)
         ENDIF

c------------Number of roots and weights equal to  14--------------------
         IF(NGQP.eq.14) then
         call root14(TVAL, UVAL, RTS, WTS, NT)
         ENDIF


c------------------------------------------------------------------------
c            DO I=1,NTGQP

c            write(6,*) RTS(I), I
c            ENDDO

C
C             ...proceed according to the number of roots and weights
C                required.


c 6000    DO N = 1,NT
c            T = TVAL (N)
c            write(6,*) 'U=',UVAL(N)
c
c         END DO

c         CALL  ERD__RYS_X_ROOTS_WEIGHTS_YUKAWA
c     +
c     +              ( NT,NTGQP,
c     +                NGQP,NMOM,
c     +                TVAL,UVAL,RYSZERO,
c     +                A,B,
c     +                MOM,
c     +                DIA,OFF,
c     +                ROW1,ROW2,
c     +
c     +                          RTS,
c     +                          WTS )
c     +
c     +
c
C
C             ...ready!
C
C
         RETURN
         END
