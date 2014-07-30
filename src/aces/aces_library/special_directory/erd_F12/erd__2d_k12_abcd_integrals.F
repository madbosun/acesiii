         SUBROUTINE  ERD__2D_K12_ABCD_INTEGRALS
     +
     +                    ( SHELLA,SHELLB,SHELLC,SHELLD,
     +                      DERA,DERB,DERC,DERD,
     +                      NGQEXQ,
     +                      EXP2A,EXP2B,EXP2C,EXP2D,
     +                      DER1,DER2,DER3,DER4,
     +                      INT2D,NNN,SWAP12,SWAP34,
     +
     +                               OUT2D,OUT2D_2 )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_DERV_ABCD_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                the input 2D ABCD integrals:
C
C         I'(n,a,b,c,d) = delta (DERA,1) *
C                         [-a*I(n,a-1,b,c,d) + 2*expa(n)*I(n,a+1,b,c,d)]
C                       + delta (DERB,1) *
C                         [-b*I(n,a,b-1,c,d) + 2*expb(n)*I(n,a,b+1,c,d)]
C                       + delta (DERC,1) *
C                         [-c*I(n,a,b,c-1,d) + 2*expc(n)*I(n,a,b,c+1,d)]
C                       + delta (DERD,1) *
C                         [-d*I(n,a,b,c,d-1) + 2*expd(n)*I(n,a,b,c,d+1)]
C
C                and returns the result in a separate array.
C
C                The derivatives of the 2D integrals are calculated for
C                all roots and the present set of exponent quadruplets.
C                The values of DERA,DERB,DERC and DERD can be only 1
C                or 0. A set of 1's among these variables indicates
C                that these centers were identified as equal.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for centers
C                                   x = A,B,C,D after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   centers x = A,B,C,D. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    EXP2x       =  the NGQEXQ exponents x 2 for all
C                                   centers x = A,B,C,D in the
C                                   appropriate order.
C                    INT2D       =  all input 2D ABCD integrals before
C                                   differentiation.
C
C
C                  Output:
C
C                    OUT2D       =  all differentiated 2D ABCD
C                                   integrals
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
         IMPLICIT  NONE

         LOGICAL   FIRST
         LOGICAL   SWAP12,SWAP34

         INTEGER   A,B,C,D,N,NNN
         INTEGER   AM,AP,BM,BP,CM,CP,DM,DP
         INTEGER   DERA,DERB,DERC,DERD
         INTEGER   NGQEXQ,DER1,DER2,DER3,DER4
         INTEGER   SHELLA,SHELLB,SHELLC,SHELLD

         DOUBLE PRECISION  F, FUCK1,FUCK2,FUCK3,FUCK4
         DOUBLE PRECISION  NEGONE,NEGTWO

         DOUBLE PRECISION  EXP2A (1:NGQEXQ)
         DOUBLE PRECISION  EXP2B (1:NGQEXQ)
         DOUBLE PRECISION  EXP2C (1:NGQEXQ)
         DOUBLE PRECISION  EXP2D (1:NGQEXQ)

         DOUBLE PRECISION  INT2D (1:NGQEXQ,0:SHELLA+DERA,
     +                                     0:SHELLB+DERB,
     +                                     0:SHELLC+DERC,
     +                                     0:SHELLD+DERD)
         DOUBLE PRECISION  OUT2D (1:NGQEXQ,0:SHELLA,
     +                                     0:SHELLB,
     +                                     0:SHELLC,
     +                                     0:SHELLD)

         DOUBLE PRECISION  OUT2D_2 (1:NGQEXQ,0:SHELLA,
     +                                       0:SHELLB,
     +                                       0:SHELLC,
     +                                       0:SHELLD)
      



         PARAMETER  (NEGONE = -1.D0)
         PARAMETER  (NEGTWO = -2.D0)
C-------------------------------------------------------------------------
C  Since we calculate integrals of [T12,f12] first derivatives of        | 
C  C and D centers should be taken with minus sign. Second derivatives   |
C  will be generated as a derivatives from first derivatives. So we need |
C  minus sign only once. Variable FUCK will take care of that minus sign |
C-------------------------------------------------------------------------

         FUCK1 = 1.d0
         FUCK2 = 1.d0
         FUCK3 = 1.d0
         FUCK4 = 1.d0 
           
         IF(NNN.eq.1) then 

               IF(.not.SWAP12) THEN
               FUCK1=   1.d0    
               FUCK2=  -1.d0
                 ELSE
               FUCK1=  -1.d0
               FUCK2=   1.d0
               ENDIF

               IF(.not.SWAP34) THEN
               FUCK3 =  1.d0
               FUCK4=  -1.d0
                 ELSE
               FUCK3=  -1.d0
               FUCK4=   1.d0
               ENDIF


         ENDIF
C------------------------------------------------------------------------
C
C
C             ...derivative on center A (if any).
C
         

c         write(6,*) '==================================='
c         write(6,*) DERA,DERB,DERC,DERD
c         write(6,*) '==================================='
C
c    

         FIRST = .TRUE.

c--------Jump to the corresponding derivative-------------------------
         IF(DER1.eq.1) goto 1111
         IF(DER2.eq.1) goto 6666
         IF(DER3.eq.1) goto 7777
         IF(DER4.eq.1) goto 8888
c---------------------------------------------------------------------

1111     continue   
         IF (DERA.EQ.1) THEN

           IF (SHELLA.EQ.0) THEN

             DO 100 D = 0,SHELLD
             DO 100 C = 0,SHELLC
             DO 100 B = 0,SHELLB
             DO 100 N = 1,NGQEXQ
                OUT2D (N,0,B,C,D) = FUCK1*EXP2A (N) * INT2D (N,1,B,C,D)
  100        CONTINUE

           ELSE IF (SHELLA.EQ.1) THEN
             DO 110 D = 0,SHELLD
             DO 110 C = 0,SHELLC
             DO 110 B = 0,SHELLB
             DO 110 N = 1,NGQEXQ
                OUT2D (N,0,B,C,D) = FUCK1*EXP2A (N) * INT2D (N,1,B,C,D)
                OUT2D (N,1,B,C,D) = FUCK1*(EXP2A (N) * INT2D (N,2,B,C,D)
     +                                        - INT2D (N,0,B,C,D) )
  110        CONTINUE

           ELSE

             DO 120 D = 0,SHELLD
             DO 120 C = 0,SHELLC
             DO 120 B = 0,SHELLB
                DO 122 N = 1,NGQEXQ
                 OUT2D (N,0,B,C,D) = FUCK1*EXP2A (N) * INT2D (N,1,B,C,D)
                 OUT2D (N,1,B,C,D) = FUCK1*(EXP2A(N) * INT2D (N,2,B,C,D)
     +                                           - INT2D (N,0,B,C,D))
  122           CONTINUE
                F = NEGTWO
                DO 124 A = 2,SHELLA
                   AM = A - 1
                   AP = A + 1
                   DO 126 N = 1,NGQEXQ
               OUT2D (N,A,B,C,D) = FUCK1*(EXP2A (N) * INT2D (N,AP,B,C,D)
     +                                        + F * INT2D (N,AM,B,C,D) )
  126              CONTINUE
                   F = F + NEGONE
  124           CONTINUE
  120        CONTINUE

           END IF

           FIRST = .FALSE.

         END IF

         goto 9999
C
C
C             ...derivative on center B (if any).
C
C
6666     continue
         IF (DERB.EQ.1) THEN

          IF (FIRST) THEN

           IF (SHELLB.EQ.0) THEN

             DO 200 D = 0,SHELLD
             DO 200 C = 0,SHELLC
             DO 200 A = 0,SHELLA
             DO 200 N = 1,NGQEXQ
                OUT2D (N,A,0,C,D) = FUCK2*EXP2B (N) * INT2D (N,A,1,C,D)
  200        CONTINUE

           ELSE IF (SHELLB.EQ.1) THEN

             DO 210 D = 0,SHELLD
             DO 210 C = 0,SHELLC
                DO 212 A = 0,SHELLA
                DO 212 N = 1,NGQEXQ
                OUT2D (N,A,0,C,D) = FUCK2*EXP2B (N) * INT2D (N,A,1,C,D)
  212           CONTINUE
                DO 214 A = 0,SHELLA
                DO 214 N = 1,NGQEXQ
                OUT2D (N,A,1,C,D) = FUCK2*(EXP2B (N) * INT2D (N,A,2,C,D)
     +                                            - INT2D (N,A,0,C,D) )
  214           CONTINUE
  210        CONTINUE

           ELSE

             DO 220 D = 0,SHELLD
             DO 220 C = 0,SHELLC
                DO 222 A = 0,SHELLA
                DO 222 N = 1,NGQEXQ
                OUT2D (N,A,0,C,D) = FUCK2*EXP2B (N) * INT2D (N,A,1,C,D)
  222           CONTINUE
                DO 224 A = 0,SHELLA
                DO 224 N = 1,NGQEXQ
                OUT2D (N,A,1,C,D) = FUCK2*(EXP2B (N) * INT2D (N,A,2,C,D)
     +                                            - INT2D (N,A,0,C,D) )
  224           CONTINUE
                F = NEGTWO
                DO 226 B = 2,SHELLB
                   BM = B - 1
                   BP = B + 1
                   DO 228 A = 0,SHELLA
                   DO 228 N = 1,NGQEXQ
               OUT2D (N,A,B,C,D) = FUCK2*(EXP2B (N) * INT2D (N,A,BP,C,D)
     +                                        + F * INT2D (N,A,BM,C,D) )
  228              CONTINUE
                   F = F + NEGONE
  226           CONTINUE
  220        CONTINUE

           END IF

           FIRST = .FALSE.

          ELSE
 
           IF (SHELLB.EQ.0) THEN

             DO 230 D = 0,SHELLD
             DO 230 C = 0,SHELLC
             DO 230 A = 0,SHELLA
             DO 230 N = 1,NGQEXQ
                OUT2D (N,A,0,C,D) = OUT2D (N,A,0,C,D) +
     +                              FUCK2*EXP2B (N) * INT2D (N,A,1,C,D)
  230        CONTINUE

           ELSE IF (SHELLB.EQ.1) THEN

             DO 240 D = 0,SHELLD
             DO 240 C = 0,SHELLC
                DO 242 A = 0,SHELLA
                DO 242 N = 1,NGQEXQ
                   OUT2D (N,A,0,C,D) = OUT2D (N,A,0,C,D) +
     +                               FUCK2*EXP2B (N) * INT2D (N,A,1,C,D)
  242           CONTINUE
                DO 244 A = 0,SHELLA
                DO 244 N = 1,NGQEXQ
                   OUT2D (N,A,1,C,D) = OUT2D (N,A,1,C,D) +
     +                              FUCK2*(EXP2B (N) * INT2D (N,A,2,C,D)
     +                                           - INT2D (N,A,0,C,D) )
  244           CONTINUE
  240        CONTINUE

           ELSE

             DO 250 D = 0,SHELLD
             DO 250 C = 0,SHELLC
                DO 252 A = 0,SHELLA
                DO 252 N = 1,NGQEXQ
                   OUT2D (N,A,0,C,D) = OUT2D (N,A,0,C,D) +
     +                               FUCK2*EXP2B (N) * INT2D (N,A,1,C,D)
  252           CONTINUE
                DO 254 A = 0,SHELLA
                DO 254 N = 1,NGQEXQ
                   OUT2D (N,A,1,C,D) = OUT2D (N,A,1,C,D) +
     +                            FUCK2*(EXP2B (N) * INT2D (N,A,2,C,D)
     +                                           - INT2D (N,A,0,C,D) )
  254           CONTINUE
                F = NEGTWO
                DO 256 B = 2,SHELLB
                   BM = B - 1
                   BP = B + 1
                   DO 258 A = 0,SHELLA
                   DO 258 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = OUT2D (N,A,B,C,D) +
     +                             FUCK2*(EXP2B (N) * INT2D (N,A,BP,C,D)
     +                                        + F * INT2D (N,A,BM,C,D) )
  258              CONTINUE
                   F = F + NEGONE
  256           CONTINUE
  250        CONTINUE

           END IF

          END IF

         END IF
         goto 9999
C
C
C             ...derivative on center C (if any).
C
C
c 
7777     continue  
         IF (DERC.EQ.1) THEN

          IF (FIRST) THEN
c             write(6,*) '@@@@@@@@@@@@@@@@@@'

           IF (SHELLC.EQ.0) THEN

             DO 300 D = 0,SHELLD
             DO 300 B = 0,SHELLB
             DO 300 A = 0,SHELLA
             DO 300 N = 1,NGQEXQ
                OUT2D (N,A,B,0,D) = FUCK3*EXP2C (N) * INT2D (N,A,B,1,D)
  300        CONTINUE

           ELSE IF (SHELLC.EQ.1) THEN
c                write(6,*) '@@@@@@@@@@@@@@@@&&&&&&&&&&&&&&&&&@@'
             DO 310 D = 0,SHELLD
                DO 312 B = 0,SHELLB
                DO 312 A = 0,SHELLA
                DO 312 N = 1,NGQEXQ
                 OUT2D (N,A,B,0,D) = FUCK3*EXP2C (N) * INT2D (N,A,B,1,D)
  312           CONTINUE
                DO 314 B = 0,SHELLB
                DO 314 A = 0,SHELLA
                DO 314 N = 1,NGQEXQ
                OUT2D (N,A,B,1,D) = FUCK3*(EXP2C (N) * INT2D (N,A,B,2,D)
     +                                           - INT2D (N,A,B,0,D) )
  314           CONTINUE
  310        CONTINUE

           ELSE

             DO 320 D = 0,SHELLD
                DO 322 B = 0,SHELLB
                DO 322 A = 0,SHELLA
                DO 322 N = 1,NGQEXQ
                 OUT2D (N,A,B,0,D) = FUCK3*EXP2C (N) * INT2D (N,A,B,1,D)
  322           CONTINUE
                DO 324 B = 0,SHELLB
                DO 324 A = 0,SHELLA
                DO 324 N = 1,NGQEXQ
                OUT2D (N,A,B,1,D) = FUCK3*(EXP2C (N) * INT2D (N,A,B,2,D)
     +                                             - INT2D (N,A,B,0,D) )
  324           CONTINUE
                F = NEGTWO
                DO 326 C = 2,SHELLC
                   CM = C - 1
                   CP = C + 1
                   DO 328 B = 0,SHELLB
                   DO 328 A = 0,SHELLA
                   DO 328 N = 1,NGQEXQ
               OUT2D (N,A,B,C,D) = FUCK3*(EXP2C (N) * INT2D (N,A,B,CP,D)
     +                                        + F * INT2D (N,A,B,CM,D) )
  328              CONTINUE
                   F = F + NEGONE
  326           CONTINUE
  320        CONTINUE

           END IF

           FIRST = .FALSE.

          ELSE

           IF (SHELLC.EQ.0) THEN

             DO 330 D = 0,SHELLD
             DO 330 B = 0,SHELLB
             DO 330 A = 0,SHELLA
             DO 330 N = 1,NGQEXQ
                OUT2D (N,A,B,0,D) = OUT2D (N,A,B,0,D) +
     +                        FUCK3* EXP2C (N) * INT2D (N,A,B,1,D)
  330        CONTINUE

           ELSE IF (SHELLC.EQ.1) THEN

             DO 340 D = 0,SHELLD
                DO 342 B = 0,SHELLB
                DO 342 A = 0,SHELLA
                DO 342 N = 1,NGQEXQ
                   OUT2D (N,A,B,0,D) = OUT2D (N,A,B,0,D) +
     +                        FUCK3* EXP2C (N) * INT2D (N,A,B,1,D)
  342           CONTINUE
                DO 344 B = 0,SHELLB
                DO 344 A = 0,SHELLA
                DO 344 N = 1,NGQEXQ
                   OUT2D (N,A,B,1,D) = OUT2D (N,A,B,1,D) +
     +                        FUCK3*(EXP2C (N) * INT2D (N,A,B,2,D)
     +                                         - INT2D (N,A,B,0,D) )
  344           CONTINUE
  340        CONTINUE

           ELSE

             DO 350 D = 0,SHELLD
                DO 352 B = 0,SHELLB
                DO 352 A = 0,SHELLA
                DO 352 N = 1,NGQEXQ
                   OUT2D (N,A,B,0,D) = OUT2D (N,A,B,0,D) +
     +                            FUCK3* EXP2C (N) * INT2D (N,A,B,1,D)
  352           CONTINUE
                DO 354 B = 0,SHELLB
                DO 354 A = 0,SHELLA
                DO 354 N = 1,NGQEXQ
                   OUT2D (N,A,B,1,D) = OUT2D (N,A,B,1,D) +
     +                             FUCK3* (EXP2C (N) * INT2D (N,A,B,2,D)
     +                                         - INT2D (N,A,B,0,D) )
  354           CONTINUE
                F = NEGTWO
                DO 356 C = 2,SHELLC
                   CM = C - 1
                   CP = C + 1
                   DO 358 B = 0,SHELLB
                   DO 358 A = 0,SHELLA
                   DO 358 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = OUT2D (N,A,B,C,D) +
     +                            FUCK3* (EXP2C (N) * INT2D (N,A,B,CP,D)
     +                                        + F * INT2D (N,A,B,CM,D) )
  358              CONTINUE
                   F = F + NEGONE
  356           CONTINUE
  350        CONTINUE

           END IF

          END IF

         END IF
         goto 9999
C
C
C             ...derivative on center D (if any).
C
C
8888     continue    
          

         IF (DERD.EQ.1) THEN


          IF (FIRST) THEN
 
           IF (SHELLD.EQ.0) THEN
             DO 400 C = 0,SHELLC
             DO 400 B = 0,SHELLB
             DO 400 A = 0,SHELLA
             DO 400 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = FUCK4*EXP2D (N) * INT2D (N,A,B,C,1)
  400        CONTINUE

           ELSE IF (SHELLD.EQ.1) THEN

             
             DO 410 C = 0,SHELLC
             DO 410 B = 0,SHELLB
             DO 410 A = 0,SHELLA
             DO 410 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = FUCK4*EXP2D (N) * INT2D (N,A,B,C,1)

  410        CONTINUE
             DO 412 C = 0,SHELLC
             DO 412 B = 0,SHELLB
             DO 412 A = 0,SHELLA
             DO 412 N = 1,NGQEXQ
           OUT2D (N,A,B,C,1) = FUCK4*(EXP2D (N) * INT2D (N,A,B,C,2)
     +       *1.d0 
     +                                       - INT2D (N,A,B,C,0) )


  412        CONTINUE

           ELSE
             DO 420 C = 0,SHELLC
             DO 420 B = 0,SHELLB
             DO 420 A = 0,SHELLA
             DO 420 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = FUCK4*EXP2D (N) * INT2D (N,A,B,C,1)
  420        CONTINUE
             DO 422 C = 0,SHELLC
             DO 422 B = 0,SHELLB
             DO 422 A = 0,SHELLA
             DO 422 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = FUCK4*(EXP2D (N) * INT2D (N,A,B,C,2)
     +                                        - INT2D (N,A,B,C,0) )
  422        CONTINUE
             F = NEGTWO
             DO 424 D = 2,SHELLD
                DM = D - 1
                DP = D + 1
                DO 426 C = 0,SHELLC
                DO 426 B = 0,SHELLB
                DO 426 A = 0,SHELLA
                DO 426 N = 1,NGQEXQ
               OUT2D (N,A,B,C,D) = FUCK4*(EXP2D (N) * INT2D (N,A,B,C,DP)
     +                                       + F * INT2D (N,A,B,C,DM) )
  426           CONTINUE
                F = F + NEGONE
  424        CONTINUE

           END IF

          ELSE

           IF (SHELLD.EQ.0) THEN

             DO 430 C = 0,SHELLC
             DO 430 B = 0,SHELLB
             DO 430 A = 0,SHELLA
             DO 430 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = OUT2D (N,A,B,C,0) +
     +                          FUCK4*EXP2D (N) * INT2D (N,A,B,C,1)
  430        CONTINUE

           ELSE IF (SHELLD.EQ.1) THEN

             DO 440 C = 0,SHELLC
             DO 440 B = 0,SHELLB
             DO 440 A = 0,SHELLA
             DO 440 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = OUT2D (N,A,B,C,0) +
     +                          FUCK4*EXP2D (N) * INT2D (N,A,B,C,1)
  440        CONTINUE
             DO 442 C = 0,SHELLC
             DO 442 B = 0,SHELLB
             DO 442 A = 0,SHELLA
             DO 442 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = OUT2D (N,A,B,C,1) +
     +                       FUCK4*(EXP2D (N) * INT2D (N,A,B,C,2)
     +                                      - INT2D (N,A,B,C,0) )
  442        CONTINUE

           ELSE

             DO 450 C = 0,SHELLC
             DO 450 B = 0,SHELLB
             DO 450 A = 0,SHELLA
             DO 450 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = OUT2D (N,A,B,C,0) +
     +                          FUCK4*EXP2D (N) * INT2D (N,A,B,C,1)
  450        CONTINUE
             DO 452 C = 0,SHELLC
             DO 452 B = 0,SHELLB
             DO 452 A = 0,SHELLA
             DO 452 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = OUT2D (N,A,B,C,1) +
     +                         FUCK4*(EXP2D (N) * INT2D (N,A,B,C,2) 
     +                                      - INT2D (N,A,B,C,0) )
  452        CONTINUE
             F = NEGTWO
             DO 454 D = 2,SHELLD
                DM = D - 1
                DP = D + 1
                DO 456 C = 0,SHELLC
                DO 456 B = 0,SHELLB
                DO 456 A = 0,SHELLA
                DO 456 N = 1,NGQEXQ
                   OUT2D (N,A,B,C,D) = OUT2D (N,A,B,C,D) +
     +                         FUCK4*( EXP2D (N) * INT2D (N,A,B,C,DP)
     +                                     + F * INT2D (N,A,B,C,DM) )
  456           CONTINUE
                F = F + NEGONE
  454        CONTINUE

           END IF

          END IF

         END IF
C
9999      continue
C
C             ...ready!
C
C--------Copy second derivatives into the scratch array---------------

         IF(NNN.eq.2) then

          IF(DER1.eq.1) then
c----------Accumulate second derivatives in OUT2D_2 array ------------
             DO D = 0,SHELLD
             DO C = 0,SHELLC
             DO B = 0,SHELLB
             DO A = 0,SHELLA 
              DO N = 1,NGQEXQ
              OUT2D_2 (N,A,B,C,D) = OUT2D (N,A,B,C,D)
              ENDDO
             ENDDO
             ENDDO
             ENDDO
             ENDDO
            
          ELSE
              
             DO D = 0,SHELLD
             DO C = 0,SHELLC
             DO B = 0,SHELLB
             DO A = 0,SHELLA
              DO N = 1,NGQEXQ
              OUT2D_2 (N,A,B,C,D) = OUT2D_2 (N,A,B,C,D) + 
     &                                     OUT2D (N,A,B,C,D)
              ENDDO
             ENDDO
             ENDDO
             ENDDO
             ENDDO
c------------Copy the sum of second derivatives to OUT2D array----------
          ENDIF

          IF(DER4.eq.1) then
             DO D = 0,SHELLD
             DO C = 0,SHELLC
             DO B = 0,SHELLB
             DO A = 0,SHELLA
              DO N = 1,NGQEXQ
              OUT2D (N,A,B,C,D) = OUT2D_2 (N,A,B,C,D)
              ENDDO
             ENDDO
             ENDDO
             ENDDO
             ENDDO
      

          ENDIF
 
         ENDIF

         RETURN
         END
