         SUBROUTINE  ERD__K12_INT2D_TO_AB00
     +
     +                    ( SHELLA,SHELLB,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZA,NXYZB,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_AB00
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris:
C
C                  [AB|00],[A0|B0],[A0|0B],[0A|B0],[0A|0B] or [00|AB]
C
C                adding up the contributions from all the respective
C                derivative 2D integrals:
C
C                        AB00,A0B0,A00B,0AB0,0A0B or 00AB
C
C                Simplified version of the general ABCD routine to
C                reduce loop overheads for those cases where there are
C                at least two s-shells. For comments and details see
C                the general ABCD routine.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x = A and B
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZx       =  # of cartesian monomials for
C                                   x = A and B shells
C                    INT2Dx      =  all current 2D AB00/A0B0/A00B/0AB0/
C                                   0A0B/00AB derivative integrals for
C                                   each cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D AB00/A0B0/A00B/0AB0/0A0B/00AB
C                                   derivative integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [AB|00]/[A0|B0]/[A0|0B]/[0A|B0]/
C                                   [0A|0B]/[00|AB] derivative
C                                   integrals corresponding to all
C                                   current exponent quadruplets
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
         IMPLICIT    NONE

         LOGICAL     DIFFY,DIFFZ

         INTEGER     I,J,M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB
         INTEGER     XA,YA,ZA,XB,YB,ZB
         INTEGER     XAP,XBP
         INTEGER     YAMAX,YBMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO,STGEXP

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ),RTS(1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ),UVAL(1:NEXQ)

         DOUBLE PRECISION  BATCH  (1:NEXQ,1:NXYZA,1:NXYZB)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLA,0:SHELLB)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to number of roots.
C
C
         GOTO  (1,2,3,4,5,6,7,8,9,10)  MIN (NGQP,10)
C
C
C                       ********************
C                       *  # of roots = 1  *
C                       ********************
C
C
    1    XBP = 0
         DO 100  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 102 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 110 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 112 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M)
                         END DO
                     ELSE
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M) * INT2DY (M,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = TEMP2 (M)*
     &                          (2.d0*UVAL(M)/STGEXP)*(1.d0-RTS(M))
                         END DO
                     ELSE
                         DO M = 1,NEXQ
                           BATCH (M,I,J) = TEMP2 (M) * INT2DZ (M,ZA,ZB)*
     &                          (2.d0*UVAL(M)/STGEXP)*(1.d0-RTS(M))
                         END DO
                     END IF

  112             CONTINUE
  110          CONTINUE
               XAP = I
  102       CONTINUE
            XBP = J
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XBP = 0
         DO 200  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 202 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 210 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 212 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                  + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 2
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 2
                         END DO
                     END IF

  212             CONTINUE
  210          CONTINUE
               XAP = I
  202       CONTINUE
            XBP = J
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XBP = 0
         DO 300  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 302 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 310 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 312 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                   + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                   + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                      )*(2.d0*UVAL(M)/STGEXP) 
                            R = R + 3
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 3
                         END DO
                     END IF

  312             CONTINUE
  310          CONTINUE
               XAP = I
  302       CONTINUE
            XBP = J
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XBP = 0
         DO 400  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 402 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 410 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 412 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                    + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                    + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                    + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 4
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA,ZB)*(1.d0-RTS(R+3))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 4
                         END DO
                     END IF

  412             CONTINUE
  410          CONTINUE
               XAP = I
  402       CONTINUE
            XBP = J
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XBP = 0
         DO 500  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 502 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 510 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 512 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                  + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                  + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                  + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                  + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 5
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA,ZB)*(1.d0-RTS(R+3))
     +                                      + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA,ZB)*(1.d0-RTS(R+4))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 5
                         END DO
                     END IF

  512             CONTINUE
  510          CONTINUE
               XAP = I
  502       CONTINUE
            XBP = J
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XBP = 0
         DO 600  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 602 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 610 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 612 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                  + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                  + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                  + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                  + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                  + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                      )*(2.d0*UVAL(M)/STGEXP) 
                            R = R + 6
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA,ZB)*(1.d0-RTS(R+3))
     +                                      + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA,ZB)*(1.d0-RTS(R+4))
     +                                      + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA,ZB)*(1.d0-RTS(R+5))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 6
                         END DO
                     END IF

  612             CONTINUE
  610          CONTINUE
               XAP = I
  602       CONTINUE
            XBP = J
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XBP = 0
         DO 700  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 702 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 710 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 712 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                  + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                  + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                  + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                  + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                  + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                  + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 7
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA,ZB)*(1.d0-RTS(R+3))
     +                                      + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA,ZB)*(1.d0-RTS(R+4))
     +                                      + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA,ZB)*(1.d0-RTS(R+5))
     +                                      + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA,ZB)*(1.d0-RTS(R+6))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 7
                         END DO
                     END IF

  712             CONTINUE
  710          CONTINUE
               XAP = I
  702       CONTINUE
            XBP = J
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XBP = 0
         DO 800  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 802 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 810 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 812 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                  + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                  + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                  + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                  + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                  + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                  + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                                  + TEMP2 (R+7)*(1.d0-RTS(R+7))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 8
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA,ZB)*(1.d0-RTS(R+3))
     +                                      + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA,ZB)*(1.d0-RTS(R+4))
     +                                      + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA,ZB)*(1.d0-RTS(R+5))
     +                                      + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA,ZB)*(1.d0-RTS(R+6))
     +                                      + TEMP2 (R+7)
     +                              * INT2DZ (R+7,ZA,ZB)*(1.d0-RTS(R+7))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 8
                         END DO
                     END IF

  812             CONTINUE
  810          CONTINUE
               XAP = I
  802       CONTINUE
            XBP = J
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XBP = 0
         DO 900  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 902 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO

               J = XBP
               DO 910 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 912 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                  + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                  + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                  + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                  + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                  + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                  + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                                  + TEMP2 (R+7)*(1.d0-RTS(R+7))
     +                                  + TEMP2 (R+8)*(1.d0-RTS(R+8))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 9
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J) = ( TEMP2 (R)
     +                                  * INT2DZ (R,ZA,ZB)*(1.d0-RTS(R))
     +                                      + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA,ZB)*(1.d0-RTS(R+1))
     +                                      + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA,ZB)*(1.d0-RTS(R+2))
     +                                      + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA,ZB)*(1.d0-RTS(R+3))
     +                                      + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA,ZB)*(1.d0-RTS(R+4))
     +                                      + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA,ZB)*(1.d0-RTS(R+5))
     +                                      + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA,ZB)*(1.d0-RTS(R+6))
     +                                      + TEMP2 (R+7)
     +                              * INT2DZ (R+7,ZA,ZB)*(1.d0-RTS(R+7))
     +                                      + TEMP2 (R+8)
     +                              * INT2DZ (R+8,ZA,ZB)*(1.d0-RTS(R+8))
     +                                      )*(2.d0*UVAL(M)/STGEXP)
                            R = R + 9
                         END DO
                     END IF

  912             CONTINUE
  910          CONTINUE
               XAP = I
  902       CONTINUE
            XBP = J
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loops over x,x-pairs. No skipping of
C                x,x-contribution of 0,0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    XBP = 0
         DO 1000  XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 1002 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB)
               END DO
C
C
C             ...inner loops over y,y-pairs. Skip the multiplication
C                of y,y-contributions, if no y-coordinate derivative
C                was formed and we have a 0,0-pair, as then the 2DY
C                integrals are equal to 1.
C
C
               J = XBP
               DO 1010 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 1012 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB)
                         END DO
                     END IF
C
C
C             ...skip multiplication of z,z-contributions, if we
C                have a 0,0-pair and no derivations were performed
C                on the z-coordinate, as then the 2DZ integrals
C                are equal to 1. All info concerning all three x,x-,
C                y,y- and z,z-pairs have been collected for all
C                exponent quadruplets at once. Sum up the 2D X,Y,Z
C                integral products to the appropriate places of the
C                batch.
C
C
                     IF (.NOT.DIFFZ .AND. ZA+ZB.EQ.0) THEN
                         R = 0
                         DO M = 1,NEXQ
                            SUM = ZERO
                            DO N = 1,NGQP
                               SUM = SUM + TEMP2 (R+N)*(1.d0-RTS(R+N))
                            END DO
                            R = R + NGQP
                            BATCH (M,I,J) = SUM*(2.d0*UVAL(M)/STGEXP)
                         END DO
                     ELSE
                         R = 0
                         DO M = 1,NEXQ
                            SUM = ZERO
                            DO N = 1,NGQP
                               SUM = SUM + TEMP2 (R+N)
     +                            * INT2DZ (R+N,ZA,ZB)*(1.d0-RTS(R+N))
                            END DO
                            R = R + NGQP
                            BATCH (M,I,J) = SUM*(2.d0*UVAL(M)/STGEXP)
                         END DO
                     END IF
C
C
C             ...next z,z-pair and y,y-pair.
C
C
 1012             CONTINUE
 1010          CONTINUE
C
C
C             ...next x,x-pair.
C
C
               XAP = I
 1002       CONTINUE
            XBP = J
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
