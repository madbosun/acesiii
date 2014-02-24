         SUBROUTINE  OED__NAI_1D_DERV_INTEGRALS
     +
     +                    ( MIJ,MGQPIJ,MGIJCEN,
     +                      NGQP,NGQPCEN,
     +                      SHELLA,SHELLB,
     +                      DERA,DERB,DERC,
     +                      EXP2A,EXP2B,EXP2AB,
     +                      IXC,P,C,
     +                      RTS,
     +                      TEMP,
     +                      INT1D,
     +
     +                               OUT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_1D_DERV_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                the input 1D nuclear attraction integrals:
C
C            I'(n,a,b) = delta (DERA,1) *
C                              [-a*I(n,a-1,b) + 2*expa(n)*I(n,a+1,b)]
C                      + delta (DERB,1) *
C                              [-b*I(n,a,b-1) + 2*expb(n)*I(n,a,b+1)]
C                      + delta (DERC,1) * r
C                              [2p(P-C)I(n,a,b) + a*I(n,a-1,b)
C                                               + b*I(n,a,b-1)]
C
C                and returns the result in a separate array.
C
C                The derivatives of the 1D integrals are calculated for
C                all roots, all nuclear centers and the present set of
C                exponent pairs. The values of DERA,DERB and DERC can
C                be only 1 or 0. A set of 1's among these variables
C                indicates that these centers were identified as equal.
C
C
C                  Input:
C
C                    MIJ         =  current # of ij primitive index
C                                   pairs corresponding to the
C                                   contracted shell pairs A,B
C                    MGQPIJ      =  # of roots times # of ij primitive
C                                   index pairs
C                    MGIJCEN     =  # of roots times # of ij primitive
C                                   index pairs times # of nuclear
C                                   attraction centers
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NGQPCEN     =  # of roots times # of nuclear
C                                   attraction centers
C                    SHELLx      =  maximum shell type for centers
C                                   x = A and B after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   centers x = A,B,C. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx
C                    EXP2x       =  the totality of all MGIJCEN
C                                   exponents x 2 for both centers
C                                   x = A and B in the appropriate
C                                   order
C                    EXP2AB      =  the MIJ distinct exponent sums
C                                   between both centers x = A and B.
C                    IXC         =  contains the index value of the
C                                   nuclear attraction center C to be
C                                   derivated. This index is necessary
C                                   for locating the appropriate
C                                   sections of the 1D integrals for
C                                   differentiation
C                    P           =  the MIJ coordinates for the
C                                   gaussian product centers P=A+B
C                    C           =  value of the coordinate for the
C                                   IXC-th nuclear attraction center
C                    RTS         =  all current MGIJCEN quadrature
C                                   roots
C                    TEMP        =  scratch array that will hold all
C                                   necessary 2pr*(P-C) values, where
C                                   r are all those MGQPIJ roots
C                                   corresponding to the IXC-th nuclear
C                                   attraction center
C                    INT1D       =  all input 1D nuclear attraction
C                                   integrals before differentiation.
C
C
C                  Output:
C
C                    OUT1D       =  all differentiated 1D nuclear
C                                   attraction integrals
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
         IMPLICIT  NONE

         LOGICAL   FIRST

         INTEGER   A,B,M,N
         INTEGER   AM,AP,BM,BP
         INTEGER   DERA,DERB,DERC
         INTEGER   IJ,NG
         INTEGER   IXC
         INTEGER   MIJ,MGQPIJ,MGIJCEN
         INTEGER   NGQP,NGQPCEN
         INTEGER   OFF
         INTEGER   SHELLA,SHELLB

         DOUBLE PRECISION  C,X
         DOUBLE PRECISION  F,F1,F2
         DOUBLE PRECISION  ZERO,ONE,TWO

         DOUBLE PRECISION  EXP2A  (1:MGIJCEN)
         DOUBLE PRECISION  EXP2B  (1:MGIJCEN)
         DOUBLE PRECISION  EXP2AB (1:MIJ)
         DOUBLE PRECISION  P      (1:MIJ)
         DOUBLE PRECISION  RTS    (1:MGIJCEN)
         DOUBLE PRECISION  TEMP   (1:MGQPIJ)

         DOUBLE PRECISION  INT1D (1:MGIJCEN,0:SHELLA+DERA,
     +                                      0:SHELLB+DERB)
         DOUBLE PRECISION  OUT1D (1:MGIJCEN,0:SHELLA,
     +                                      0:SHELLB)

         PARAMETER  (ZERO   =  0.D0)
         PARAMETER  (ONE    =  1.D0)
         PARAMETER  (TWO    =  2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...derivative on center A (if any).
C
C
         FIRST = .TRUE.

         IF (DERA.EQ.1) THEN

             IF (SHELLA.EQ.0) THEN

                 DO B = 0,SHELLB
                 DO N = 1,MGIJCEN
                    OUT1D (N,0,B) = EXP2A (N) * INT1D (N,1,B)
                 END DO
                 END DO

             ELSE IF (SHELLA.EQ.1) THEN

                 DO B = 0,SHELLB
                 DO N = 1,MGIJCEN
                    OUT1D (N,0,B) = EXP2A (N) * INT1D (N,1,B)
                    OUT1D (N,1,B) = EXP2A (N) * INT1D (N,2,B)
     +                                        - INT1D (N,0,B)
                 END DO
                 END DO

             ELSE

                 DO B = 0,SHELLB
                    DO N = 1,MGIJCEN
                       OUT1D (N,0,B) = EXP2A (N) * INT1D (N,1,B)
                       OUT1D (N,1,B) = EXP2A (N) * INT1D (N,2,B)
     +                                           - INT1D (N,0,B)
                    END DO
                    F = TWO
                    DO A = 2,SHELLA
                       AM = A - 1
                       AP = A + 1
                       DO N = 1,MGIJCEN
                          OUT1D (N,A,B) = EXP2A (N) * INT1D (N,AP,B)
     +                                          - F * INT1D (N,AM,B)
                       END DO
                       F = F + ONE
                    END DO
                 END DO

             END IF

             FIRST = .FALSE.

         END IF
C
C
C             ...derivative on center B (if any).
C
C
         IF (DERB.EQ.1) THEN

             IF (FIRST) THEN

                 DO A = 0,SHELLA
                 DO N = 1,MGIJCEN
                    OUT1D (N,A,0) = EXP2B (N) * INT1D (N,A,1)
                 END DO
                 END DO

                 IF (SHELLB.GT.0) THEN
                     DO A = 0,SHELLA
                     DO N = 1,MGIJCEN
                        OUT1D (N,A,1) = EXP2B (N) * INT1D (N,A,2)
     +                                            - INT1D (N,A,0)
                     END DO
                     END DO
                 END IF

                 IF (SHELLB.GT.1) THEN
                     F = TWO
                     DO B = 2,SHELLB
                        BM = B - 1
                        BP = B + 1
                        DO A = 0,SHELLA
                        DO N = 1,MGIJCEN
                           OUT1D (N,A,B) = EXP2B (N) * INT1D (N,A,BP)
     +                                           - F * INT1D (N,A,BM)
                        END DO
                        END DO
                        F = F + ONE
                     END DO
                 END IF

                 FIRST = .FALSE.

             ELSE
 
                 DO A = 0,SHELLA
                 DO N = 1,MGIJCEN
                    OUT1D (N,A,0) = OUT1D (N,A,0) +
     +                              EXP2B (N) * INT1D (N,A,1)
                 END DO
                 END DO

                 IF (SHELLB.GT.0) THEN
                     DO A = 0,SHELLA
                     DO N = 1,MGIJCEN
                        OUT1D (N,A,1) = OUT1D (N,A,1) +
     +                                  EXP2B (N) * INT1D (N,A,2)
     +                                            - INT1D (N,A,0)
                     END DO
                     END DO
                 END IF

                 IF (SHELLB.GT.1) THEN
                     F = TWO
                     DO B = 2,SHELLB
                        BM = B - 1
                        BP = B + 1
                        DO A = 0,SHELLA
                        DO N = 1,MGIJCEN
                           OUT1D (N,A,B) = OUT1D (N,A,B) +
     +                                     EXP2B (N) * INT1D (N,A,BP)
     +                                           - F * INT1D (N,A,BM)
                        END DO
                        END DO
                        F = F + ONE
                     END DO
                 END IF

             END IF

         END IF
C
C
C             ...derivative on center C (if any). Note, that only
C                one of the complete set of nuclear attraction centers
C                might contribute to the derivation. This center,
C                if present, is labeled by its set index IXC, and only
C                specific sections of the 1D integrals need to be
C                differentiated here.
C
C
         IF (DERC.EQ.1) THEN

             OFF = NGQP * (IXC - 1)

             M = 0
             N = OFF
             DO IJ = 1,MIJ
                X = (P (IJ) - C) * EXP2AB (IJ)
                DO NG = 1,NGQP
                   TEMP (M+NG) = RTS (N+NG) * X
                END DO
                M = M + NGQP
                N = N + NGQPCEN
             END DO

             IF (FIRST) THEN

                 DO B = 0,SHELLB
                 DO A = 0,SHELLA
                 DO N = 1,MGIJCEN
                    OUT1D (N,A,B) = ZERO
                 END DO
                 END DO
                 END DO

                 M = 0
                 N = OFF
                 DO IJ = 1,MIJ
                    DO NG = 1,NGQP
                       OUT1D (N+NG,0,0) = TEMP (M+NG) * INT1D (N+NG,0,0)
                    END DO
                    M = M + NGQP
                    N = N + NGQPCEN
                 END DO

                 IF (SHELLA.GT.0) THEN
                     F1 = ONE
                     DO A = 1,SHELLA
                        AM = A - 1
                        M = 0
                        N = OFF
                        DO IJ = 1,MIJ
                           DO NG = 1,NGQP
                              OUT1D (N+NG,A,0) =
     +                                  TEMP (M+NG) * INT1D (N+NG,A ,0)
     +                            + F1 * RTS (N+NG) * INT1D (N+NG,AM,0)
                           END DO
                           M = M + NGQP
                           N = N + NGQPCEN
                        END DO
                        F1 = F1 + ONE
                     END DO
                 END IF

                 IF (SHELLB.GT.0) THEN
                     F2 = ONE
                     DO B = 1,SHELLB
                        BM = B - 1
                        M = 0
                        N = OFF
                        DO IJ = 1,MIJ
                           DO NG = 1,NGQP
                              OUT1D (N+NG,0,B) =
     +                                  TEMP (M+NG) * INT1D (N+NG,0,B )
     +                            + F2 * RTS (N+NG) * INT1D (N+NG,0,BM)
                           END DO
                           M = M + NGQP
                           N = N + NGQPCEN
                        END DO

                        IF (SHELLA.GT.0) THEN
                            F1 = ONE
                            DO A = 1,SHELLA
                               AM = A - 1
                               M = 0
                               N = OFF
                               DO IJ = 1,MIJ
                                  DO NG = 1,NGQP
                                     OUT1D (N+NG,A,B) =
     +                                  TEMP (M+NG) * INT1D (N+NG,A ,B )
     +                            + F1 * RTS (N+NG) * INT1D (N+NG,AM,B )
     +                            + F2 * RTS (N+NG) * INT1D (N+NG,A ,BM)
                                  END DO
                                  M = M + NGQP
                                  N = N + NGQPCEN
                               END DO
                               F1 = F1 + ONE
                            END DO
                        END IF

                        F2 = F2 + ONE
                     END DO
                 END IF

             ELSE

                 M = 0
                 N = OFF
                 DO IJ = 1,MIJ
                    DO NG = 1,NGQP
                       OUT1D (N+NG,0,0) = OUT1D (N+NG,0,0)
     +                                  + TEMP (M+NG) * INT1D (N+NG,0,0)
                    END DO
                    M = M + NGQP
                    N = N + NGQPCEN
                 END DO

                 IF (SHELLA.GT.0) THEN
                     F1 = ONE
                     DO A = 1,SHELLA
                        AM = A - 1
                        M = 0
                        N = OFF
                        DO IJ = 1,MIJ
                           DO NG = 1,NGQP
                              OUT1D (N+NG,A,0) = OUT1D (N+NG,A,0)
     +                            +     TEMP (M+NG) * INT1D (N+NG,A ,0)
     +                            + F1 * RTS (N+NG) * INT1D (N+NG,AM,0)
                           END DO
                           M = M + NGQP
                           N = N + NGQPCEN
                        END DO
                        F1 = F1 + ONE
                     END DO
                 END IF

                 IF (SHELLB.GT.0) THEN
                     F2 = ONE
                     DO B = 1,SHELLB
                        BM = B - 1
                        M = 0
                        N = OFF
                        DO IJ = 1,MIJ
                           DO NG = 1,NGQP
                              OUT1D (N+NG,0,B) = OUT1D (N+NG,0,B)
     +                            +     TEMP (M+NG) * INT1D (N+NG,0,B )
     +                            + F2 * RTS (N+NG) * INT1D (N+NG,0,BM)
                           END DO
                           M = M + NGQP
                           N = N + NGQPCEN
                        END DO

                        IF (SHELLA.GT.0) THEN
                            F1 = ONE
                            DO A = 1,SHELLA
                               AM = A - 1
                               M = 0
                               N = OFF
                               DO IJ = 1,MIJ
                                  DO NG = 1,NGQP
                                     OUT1D (N+NG,A,B) = OUT1D (N+NG,A,B)
     +                            +     TEMP (M+NG) * INT1D (N+NG,A ,B )
     +                            + F1 * RTS (N+NG) * INT1D (N+NG,AM,B )
     +                            + F2 * RTS (N+NG) * INT1D (N+NG,A ,BM)
                                  END DO
                                  M = M + NGQP
                                  N = N + NGQPCEN
                               END DO
                               F1 = F1 + ONE
                            END DO
                        END IF

                        F2 = F2 + ONE
                     END DO
                 END IF

             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
