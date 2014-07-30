
      SUBROUTINE FOLOWTANGENT(HESMOD, DIAGHES, QSTLST_TANGENT,
     &                        SCRATCH, NX, NOPT, IMODE)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      DOUBLE PRECISION OVRLP

      PARAMETER (OVRLP_MIN = 8.0D-1)
      DATA ZERO /0.0D0/
      DATA TWO  /2.0D0/

      DIMENSION HESMOD(NOPT, NOPT), DIAGHES(NOPT, NOPT),
     &          SCRATCH(NX*NX),  QSTLST_TANGENT(NOPT)
C
C For climbing steps, Let's get the eigenvector of the Hessian that 
C has the largest overlap with the LST or QST tangent. Since the 
C tangent is approx. eigenfunction (by construction) we are basically
C following the tangent.
C
      ZOVRLP1 = ZERO
      ZOVRLP2 = ZERO

      DO 10 IOPT = 1, NOPT
         OVRLP = XDDOT(NOPT,QSTLST_TANGENT,1,DIAGHES(1,IOPT),1)
         IF (DABS(OVRLP).GT.ZOVRLP1) THEN
            ZOVRLP2 = ZOVRLP1
            ZOVRLP1 = DABS(OVRLP)
            MEFOLOW = IOPT
         ELSE IF (DABS(OVRLP) .GT. ZOVRLP2) THEN
            ZOVRLP2 = DABS(OVRLP)
         END IF
 10   CONTINUE

      IMODE = MEFOLOW

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,F12.6,1x,F12.6,1x,I3)") "The maximum overlaps = ", 
     &                                        ZOVRLP1,ZOVRLP2, MEFOLOW
      Write(6,*) 
      Write(6,"(a,I3,a,a)")"The LST/QST mode ",Imode," is being used",
     &                     " for climbing."
#endif

      RETURN
      END

