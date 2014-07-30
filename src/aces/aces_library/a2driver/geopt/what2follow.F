
      SUBROUTINE WHAT2FOLLOW(HESMOD, DIAGHES, QSTLST_TANGENT,
     &                       SCRATCH, NX, NOPT, IMODE)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      DOUBLE PRECISION OVRLP

      PARAMETER (OVRLP_MIN = 8.0D-1)
      DATA ZERO /0.0D0/
      DATA TWO  /2.0D0/

      DIMENSION HESMOD(NOPT, NOPT), DIAGHES(NOPT, NOPT),
     &          SCRATCH(NX*NX),  QSTLST_TANGENT(NOPT)
C
C Let's get the eigenvector of the Hessian that has the largest overlap
C with the LST or QST tangent. If there is no eigenvector that has
C an overlap of larger than 0.8 or twice as big as the next
C largest overlap, then follow the eigenvector corresponding to
C the smallest eigenvalue. If the user wants to follow something
C other than said eigenvectors, then there is no point of doing
C LST or QST.
C
      ZOVRLP1 = ZERO
      ZOVRLP2 = ZERO

      DO 10 IOPT = 1, NOPT
         OVRLP = DDOT(NOPT,QSTLST_TANGENT,1,DIAGHES(1,IOPT),1)
         IF (DABS(OVRLP).GT.ZOVRLP1) THEN
            ZOVRLP2 = ZOVRLP1
            ZOVRLP1 = DABS(OVRLP)
            MEFOLOW = IOPT
         ELSE IF (DABS(OVRLP) .GT. ZOVRLP2) THEN
            ZOVRLP2 = DABS(OVRLP)
         END IF
 10   CONTINUE

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,F12.6,1x,F12.6,1x,I3)") "The maximum overlaps = ", 
     &                                      ZOVRLP1,ZOVRLP2,
     &                                      MEFOLOW
#endif

      IF (ZOVRLP1 .GT. OVRLP_MIN .OR. ZOVRLP1 .GT. TWO*ZOVRLP2) THEN
         IMODE = MEFOLOW

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,I3)")"The LST/QST mode that is being followed  = ", 
     &                    Imode
#endif
C
C The pure eigenvector following logic. ie. follow the eiegenvector corresponds
C to the lowest eigenvalue (or user's choice). Find out the eigenvector
C for the lowest eigenvalue (IVEC=1). That is the default.
C
      ELSE

          DO 30 J = 1, NOPT
             SCRATCH(J) = HESMOD(J, J)
 30       CONTINUE
          Z = 100.D0
          DO 50 I = 1, NOPT
             IF (SCRATCH(I) .LT. Z) THEN
                IMODE = I
                Z = SCRATCH(I)
             END IF
 50       CONTINUE

#ifdef _DEBUG_LVL0
      Write(6,"(a,a,1x,I3)")"The LST/QST, but lowest eig. value is",
     &                   " being followed  = ", Imode
#endif

      ENDIF

      RETURN
      END

