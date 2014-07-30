         SUBROUTINE  ERD__K12_INT2D_TO_ABCD
     +
     +                    ( SHELLA,SHELLB,SHELLC,SHELLD,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZA,NXYZB,NXYZC,NXYZD,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,RTS,UVAL,STGEXP,  
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_ABCD
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris [AB|CD], adding up the contributions
C                from all the derivative 2D ABCD integrals.
C
C                The routine uses the reduced Rys multiplication scheme
C                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889.
C                This scheme reuses intermediate products between
C                2DX and 2DY integrals.
C
C                Due to the very computational intensive steps inside
C                the x,y,z loops, special sections of identical x,y,z
C                loop structures have been given for each # of roots
C                =< 9, thus saving considerable computing time over the
C                general case.
C
C                For comments on how the x,y,z loop structures are
C                coded please refer to the general root case.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x = A,B,C,D
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZx       =  # of cartesian monomials for
C                                   x = A,B,C,D shells
C                    INT2Dx      =  all current 2D ABCD derivative
C                                   integrals for each cartesian
C                                   component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D ABCD derivative integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [AB|CD] derivative integrals
C                                   corresponding to all current
C                                   exponent quadruplets
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

         INTEGER     I,J,K,L,M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER     SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER     XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
         INTEGER     XAP,XBP,XCP,XDP
         INTEGER     YAMAX,YBMAX,YCMAX,YDMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO,STGEXP,QQQ

         DOUBLE PRECISION  SCALE (1:NGQEXQ),RTS(1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ),UVAL(1:NEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH  (1:NEXQ,
     +                             1:NXYZA,1:NXYZB,1:NXYZC,1:NXYZD)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,
     +                             0:SHELLA,0:SHELLB,0:SHELLC,0:SHELLD)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,
     +                             0:SHELLA,0:SHELLB,0:SHELLC,0:SHELLD)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,
     +                             0:SHELLA,0:SHELLB,0:SHELLC,0:SHELLD)

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
    1    XDP = 0
         DO 100 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 102 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 104  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 106 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 110 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 112 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 114 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 116 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO M = 1,NEXQ
                                TEMP2 (M) = TEMP1 (M)
                             END DO
                         ELSE
                             DO M = 1,NEXQ
                                TEMP2 (M) = TEMP1 (M)
     +                                    * INT2DY (M,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = TEMP2 (M)*
     +                           (2.d0*UVAL(M)/STGEXP)*(1.d0-RTS(M))
                             END DO
                         ELSE
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) =
     +                                           TEMP2 (M)
     +                                         * INT2DZ (M,ZA,ZB,ZC,ZD)*
     +                           (2.d0*UVAL(M)/STGEXP)*(1.d0-RTS(M))
                             END DO
                         END IF

  116                  CONTINUE
  114                CONTINUE
  112              CONTINUE
  110            CONTINUE
                 XAP = I
  106          CONTINUE
               XBP = J
  104        CONTINUE
             XCP = K
  102      CONTINUE
           XDP = L
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XDP = 0
         DO 200 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 202 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 204  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 206 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 210 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 212 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 214 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 216 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                          )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 2
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                        )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 2
                             END DO
                         END IF

  216                  CONTINUE
  214                CONTINUE
  212              CONTINUE
  210            CONTINUE
                 XAP = I
  206          CONTINUE
               XBP = J
  204        CONTINUE
             XCP = K
  202      CONTINUE
           XDP = L
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XDP = 0
         DO 300 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 302 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 304  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 306 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 310 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 312 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 314 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 316 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                         )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 3
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                           )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 3
                             END DO
                         END IF

  316                  CONTINUE
  314                CONTINUE
  312              CONTINUE
  310            CONTINUE
                 XAP = I
  306          CONTINUE
               XBP = J
  304        CONTINUE
             XCP = K
  302      CONTINUE
           XDP = L
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XDP = 0
         DO 400 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 402 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 404  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 406 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 410 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 412 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 414 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 416 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                     + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                           )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 4
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = ( 
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                       + TEMP2 (R+3)
     +                        * INT2DZ (R+3,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+3))
     +                                         )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 4
                             END DO
                         END IF

  416                  CONTINUE
  414                CONTINUE
  412              CONTINUE
  410            CONTINUE
                 XAP = I
  406          CONTINUE
               XBP = J
  404        CONTINUE
             XCP = K
  402      CONTINUE
           XDP = L
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XDP = 0
         DO 500 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 502 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 504  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 506 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 510 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 512 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 514 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 516 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                     + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                     + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                         )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 5
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                       + TEMP2 (R+3)
     +                        * INT2DZ (R+3,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+3))
     +                                       + TEMP2 (R+4)
     +                        * INT2DZ (R+4,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+4))
     +                                       )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 5
                             END DO
                         END IF

  516                  CONTINUE
  514                CONTINUE
  512              CONTINUE
  510            CONTINUE
                 XAP = I
  506          CONTINUE
               XBP = J
  504        CONTINUE
             XCP = K
  502      CONTINUE
           XDP = L
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XDP = 0
         DO 600 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 602 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 604  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 606 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 610 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 612 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 614 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 616 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                     + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                     + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                     + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                          )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 6
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                       + TEMP2 (R+3)
     +                        * INT2DZ (R+3,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+3))
     +                                       + TEMP2 (R+4)
     +                        * INT2DZ (R+4,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+4))
     +                                       + TEMP2 (R+5)
     +                        * INT2DZ (R+5,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+5))
     +                                          )*(2.d0*UVAL(M)/STGEXP) 
                                R = R + 6
                             END DO
                         END IF

  616                  CONTINUE
  614                CONTINUE
  612              CONTINUE
  610            CONTINUE
                 XAP = I
  606          CONTINUE
               XBP = J
  604        CONTINUE
             XCP = K
  602      CONTINUE
           XDP = L
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XDP = 0
         DO 700 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 702 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 704  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 706 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 710 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 712 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 714 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 716 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                     + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                     + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                     + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                     + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                                         )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 7
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                       + TEMP2 (R+3)
     +                        * INT2DZ (R+3,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+3))
     +                                       + TEMP2 (R+4)
     +                        * INT2DZ (R+4,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+4))
     +                                       + TEMP2 (R+5)
     +                        * INT2DZ (R+5,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+5))
     +                                       + TEMP2 (R+6)
     +                        * INT2DZ (R+6,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+6))
     +                                         )*(2.d0*UVAL(M)/STGEXP)  
                                R = R + 7
                             END DO
                         END IF

  716                  CONTINUE
  714                CONTINUE
  712              CONTINUE
  710            CONTINUE
                 XAP = I
  706          CONTINUE
               XBP = J
  704        CONTINUE
             XCP = K
  702      CONTINUE
           XDP = L
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XDP = 0
         DO 800 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 802 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 804  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 806 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 810 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 812 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 814 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 816 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                     + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                     + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                     + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                     + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                                     + TEMP2 (R+7)*(1.d0-RTS(R+7))
     +                                       )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 8
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                       + TEMP2 (R+3)
     +                        * INT2DZ (R+3,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+3))
     +                                       + TEMP2 (R+4)
     +                        * INT2DZ (R+4,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+4))
     +                                       + TEMP2 (R+5)
     +                        * INT2DZ (R+5,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+5))
     +                                       + TEMP2 (R+6)
     +                        * INT2DZ (R+6,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+6))
     +                                       + TEMP2 (R+7)
     +                        * INT2DZ (R+7,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+7))
     +                                        )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 8
                             END DO
                         END IF

  816                  CONTINUE
  814                CONTINUE
  812              CONTINUE
  810            CONTINUE
                 XAP = I
  806          CONTINUE
               XBP = J
  804        CONTINUE
             XCP = K
  802      CONTINUE
           XDP = L
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XDP = 0
         DO 900 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 902 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 904  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 906 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO

                 L = XDP
                 DO 910 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 912 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 914 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 916 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF

                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 1
                             DO M = 1,NEXQ
                           BATCH (M,I,J,K,L) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                                     + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                                     + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                                     + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                                     + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                                     + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                                     + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                                     + TEMP2 (R+7)*(1.d0-RTS(R+7))
     +                                     + TEMP2 (R+8)*(1.d0-RTS(R+8))
     +                                          )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 9
                             END DO
                         ELSE
                             R = 1
                             DO M = 1,NEXQ
                                BATCH (M,I,J,K,L) = (
     +                                         TEMP2 (R)
     +                        * INT2DZ (R,ZA,ZB,ZC,ZD)*(1.d0-RTS(R))
     +                                       + TEMP2 (R+1)
     +                        * INT2DZ (R+1,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+1))
     +                                       + TEMP2 (R+2)
     +                        * INT2DZ (R+2,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+2))
     +                                       + TEMP2 (R+3)
     +                        * INT2DZ (R+3,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+3))
     +                                       + TEMP2 (R+4)
     +                        * INT2DZ (R+4,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+4))
     +                                       + TEMP2 (R+5)
     +                        * INT2DZ (R+5,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+5))
     +                                       + TEMP2 (R+6)
     +                        * INT2DZ (R+6,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+6))
     +                                       + TEMP2 (R+7)
     +                        * INT2DZ (R+7,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+7))
     +                                       + TEMP2 (R+8)
     +                        * INT2DZ (R+8,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+8))
     +                                           )*(2.d0*UVAL(M)/STGEXP)
                                R = R + 9
                             END DO
                         END IF

  916                  CONTINUE
  914                CONTINUE
  912              CONTINUE
  910            CONTINUE
                 XAP = I
  906          CONTINUE
               XBP = J
  904        CONTINUE
             XCP = K
  902      CONTINUE
           XDP = L
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loops over x,x,x,x-quadruples. No skipping of
C                x,x,x,x-contribution of 0,0,0,0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    XDP = 0
         DO 1000 XD = SHELLD,0,-1
           YDMAX = SHELLD - XD
           XCP = 0
           DO 1002 XC = SHELLC,0,-1
             YCMAX = SHELLC - XC
             XBP = 0
             DO 1004  XB = SHELLB,0,-1
               YBMAX = SHELLB - XB
               XAP = 0
               DO 1006 XA = SHELLA,0,-1
                 YAMAX = SHELLA - XA

                 DO M = 1,NGQEXQ
                    TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC,XD)
                 END DO
C
C
C             ...inner loops over y,y,y,y-quadruples. Skip the
C                multiplication of y,y,y,y-contributions, if no
C                y-coordinate derivative was formed and we have a
C                0,0,0,0-quadruple, as then the 2DY ABCD integrals
C                are equal to 1.
C
C
                 L = XDP
                 DO 1010 YD = YDMAX,0,-1
                   L = L + 1
                   ZD = YDMAX - YD
                   K = XCP
                   DO 1012 YC = YCMAX,0,-1
                     K = K + 1
                     ZC = YCMAX - YC
                     J = XBP
                     DO 1014 YB = YBMAX,0,-1
                       J = J + 1
                       ZB = YBMAX - YB
                       I = XAP
                       DO 1016 YA = YAMAX,0,-1
                         I = I + 1
                         ZA = YAMAX - YA

                         IF (.NOT.DIFFY .AND. YA+YB+YC+YD.EQ.0) THEN
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
                             END DO
                         ELSE
                             DO N = 1,NGQEXQ
                                TEMP2 (N) = TEMP1 (N)
     +                                    * INT2DY (N,YA,YB,YC,YD)
                             END DO
                         END IF
C
C
C             ...skip multiplication of z,z,z,z-contributions, if we
C                have a 0,0,0,0-quadruple and no derivations were
C                performed on the z-coordinate, as then the 2DZ ABCD
C                integrals are equal to 1. All info concerning all
C                three x,x,x,x-, y,y,y,y- and z,z,z,z-quadruples have
C                been collected for all exponent quadruplets at once.
C                Sum up the 2D X,Y,Z integral products to the
C                appropriate places of the [AB|CD] batch.
C
C
                         IF (.NOT.DIFFZ .AND. ZA+ZB+ZC+ZD.EQ.0) THEN
                             R = 0
                             DO M = 1,NEXQ
                                SUM = ZERO
                                DO N = 1,NGQP
                                 QQQ = TEMP2 (R+N)*(1.d0-RTS(R+N))
                                 SUM = SUM + QQQ
                                END DO
                                R = R + NGQP
                           BATCH (M,I,J,K,L) = SUM*(2.d0*UVAL(M)/STGEXP)
                             END DO
                         ELSE
                             R = 0
                             DO M = 1,NEXQ
                                SUM = ZERO
                                DO N = 1,NGQP
             QQQ = TEMP2 (R+N)*INT2DZ(R+N,ZA,ZB,ZC,ZD)*(1.d0-RTS(R+N))
                                   SUM = SUM + QQQ
                                END DO
                                R = R + NGQP
                           BATCH (M,I,J,K,L) = SUM*(2.d0*UVAL(M)/STGEXP)
                             END DO
                         END IF
C
C
C             ...next z,z,z,z-quadruple and y,y,y,y-quadruple.
C
C
 1016                  CONTINUE
 1014                CONTINUE
 1012              CONTINUE
 1010            CONTINUE
C
C
C             ...next x,x,x,x-quadruple.
C
C
                 XAP = I
 1006          CONTINUE
               XBP = J
 1004        CONTINUE
             XCP = K
 1002      CONTINUE
           XDP = L
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
