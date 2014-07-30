         SUBROUTINE  ERD__K12_INT2D_TO_A000
     +
     +                    ( SHELLA,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZA,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_A000
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris:
C
C                     [A0|00] , [0A|00] , [00|A0] or [00|0A]
C
C                adding up the contributions from all the respective
C                derivative 2D integrals:
C
C                           A000 , 0A00 , 00A0 or 000A
C
C                Simplified version of the general ABCD routine to
C                reduce loop overheads for those cases where there are
C                at least three s-shells. For comments and details see
C                the general ABCD routine.
C
C
C                  Input:
C
C                    SHELLA      =  shell type for contraction shell A
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZA       =  # of cartesian monomials for shell A
C                    INT2Dx      =  all current 2D A000/0A00/00A0/000A
C                                   derivative integrals for each
C                                   cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D A000/0A00/00A0/000A derivative
C                                   integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [A0|00]/[0A|00]/[00|A0]/[00|0A]
C                                   derivative integrals corresponding
C                                   to all current exponent quadruplets
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

         INTEGER     I,M,N,R,i1,i2
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZA
         INTEGER     SHELLA
         INTEGER     XA,YA,ZA
         INTEGER     YAMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO,STGEXP

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ),RTS(1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ),UVAL(1:NEXQ)

         DOUBLE PRECISION  BATCH  (1:NEXQ,1:NXYZA)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLA)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLA)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLA)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to number of roots.
C
c           write(6,*) STGEXP,'$$$$$$$$$$$$$'    
              
c            do i1=1,NEXQ
c            write(6,*) UVAL(i1), '!!!!!!!!',i1 
c            enddo

c             do i1=1,NEXQ
c              do i2=NXYZA
c              BATCH(i1,i2) = 0.d0
c              enddo
c             enddo 


c           do i1=1,NGQEXQ
c            write(6,*) SCALE (i1), RTS(i1),'+*******+',i1
c           enddo

c           do i1=1,NGQEXQ
c            do i2=1,SHELLA
c            write(6,*) INT2DX(i1,i2),'----',i1,i2
c            write(6,*) INT2DY(i1,i2),'<---',i1,i2
c            write(6,*) INT2DZ(i1,i2),'<<--',i1,i2
c            enddo
c            enddo

c            write(6,*) NEXQ, NXYZA,'~~~~~~~~~~'
C
         GOTO  (1,2,3,4,5,6,7,8,9,10)  MIN (NGQP,10)
C
C
C                       ********************
C                       *  # of roots = 1  *
C                       ********************
C
C
    1    I = 0
         DO 100 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 110 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO M = 1,NEXQ
                      TEMP2 (M) = TEMP1 (M)
                   END DO
               ELSE
                   DO M = 1,NEXQ
                      TEMP2 (M) = TEMP1 (M) * INT2DY (M,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   DO M = 1,NEXQ
                      BATCH (M,I) = TEMP2 (M)*
     &                  (2.d0*UVAL(M)/STGEXP)*(1.d0-RTS(M))
c                     write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   DO M = 1,NEXQ
                      BATCH (M,I) = TEMP2 (M) * INT2DZ (M,ZA)*
     &                  (2.d0*UVAL(M)/STGEXP)*(1.d0-RTS(M))
c                   write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  110       CONTINUE
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    I = 0
         DO 200 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 210 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = (  TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1)) 
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 2
c                     write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 2
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  210       CONTINUE
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    I = 0
         DO 300 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 310 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 3
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 3
c                    write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  310       CONTINUE
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    I = 0
         DO 400 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 410 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                              )*(2.d0*UVAL(M)/STGEXP) 
                      R = R + 4
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)*(1.d0-RTS(R+3))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 4
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  410       CONTINUE
  400    CONTINUE




         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    I = 0
         DO 500 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 510 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 5
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)*(1.d0-RTS(R+4))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 5
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  510       CONTINUE
  500    CONTINUE



         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    I = 0
         DO 600 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 610 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 6
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)*(1.d0-RTS(R+5))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 6
c                      write(6,*) BATCH (M,I),M,I,NGQP                  
                   END DO
               END IF

  610       CONTINUE
  600    CONTINUE




         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    I = 0
         DO 700 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 710 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                              + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 7
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)*(1.d0-RTS(R+5))
     +                              + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA)*(1.d0-RTS(R+6))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 7
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  710       CONTINUE
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    I = 0
         DO 800 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 810 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                              + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                              + TEMP2 (R+7)*(1.d0-RTS(R+7))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 8
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)*(1.d0-RTS(R+5))
     +                              + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA)*(1.d0-RTS(R+6))
     +                              + TEMP2 (R+7)
     +                              * INT2DZ (R+7,ZA)*(1.d0-RTS(R+7))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 8
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  810       CONTINUE
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    I = 0
         DO 900 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 910 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)*(1.d0-RTS(R+5))
     +                              + TEMP2 (R+6)*(1.d0-RTS(R+6))
     +                              + TEMP2 (R+7)*(1.d0-RTS(R+7))
     +                              + TEMP2 (R+8)*(1.d0-RTS(R+8))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 9
c                       write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) = ( TEMP2 (R)
     +                              * INT2DZ (R,ZA)*(1.d0-RTS(R))
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)*(1.d0-RTS(R+1))
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)*(1.d0-RTS(R+2))
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)*(1.d0-RTS(R+3))
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)*(1.d0-RTS(R+4))
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)*(1.d0-RTS(R+5))
     +                              + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA)*(1.d0-RTS(R+6))
     +                              + TEMP2 (R+7)
     +                              * INT2DZ (R+7,ZA)*(1.d0-RTS(R+7))
     +                              + TEMP2 (R+8)
     +                              * INT2DZ (R+8,ZA)*(1.d0-RTS(R+8))
     +                              )*(2.d0*UVAL(M)/STGEXP)
                      R = R + 9
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF

  910       CONTINUE
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loop over x-contribution. No skipping of
C                x-contribution of 0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    I = 0
         DO 1000 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO
C
C
C             ...inner loop over y-contribution. Skip the multiplication
C                of y-contributions, if no y-coordinate derivative
C                was formed and we have a 0-contribution, as then the
C                2DY integrals are equal to 1.
C
C
            DO 1010 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF
C
C
C             ...skip multiplication of z-contributions, if we
C                have a 0-type and no derivations were performed
C                on the z-coordinate, as then the 2DZ integrals
C                are equal to 1. All info concerning all three x-,
C                y- and z-contributions have been collected for all
C                exponent quadruplets at once. Sum up the 2D X,Y,Z
C                integral products to the appropriate places of the
C                batch.
C
C
               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 0
                   DO M = 1,NEXQ
                      SUM = ZERO
                      DO N = 1,NGQP
                         SUM = SUM + TEMP2 (R+N)*(1.d0-RTS(R+N))
                      END DO
                      R = R + NGQP
                      BATCH (M,I) = SUM*(2.d0*UVAL(M)/STGEXP)
c                      write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               ELSE
                   R = 0
                   DO M = 1,NEXQ
                      SUM = ZERO
                      DO N = 1,NGQP
                         SUM = SUM + TEMP2 (R+N) * INT2DZ (R+N,ZA)*
     &                                       (1.d0-RTS(R+N))
                      END DO
                      R = R + NGQP
                      BATCH (M,I) = SUM*(2.d0*UVAL(M)/STGEXP)
                      
c                     write(6,*) BATCH (M,I),M,I,NGQP
                   END DO
               END IF
C
C
C             ...next z- and y-contributions.
C
C
 1010       CONTINUE
C
C
C             ...next x-contribution.
C
C
 1000    CONTINUE
C
C
C             ...ready!
C

11111    continue

         RETURN
         END
