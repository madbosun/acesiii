
C FORMS THE QUANTITY:
C
C         E = SUM  [E(1) + TAU(i) * [E(i)-E(1)]]   
C
C FOR QCISD, LCCD AND CCD, ONLY THE T2 PART OF THE ENERGY IS ADDED IN,
C WHILE THE T1**2 PIECE IS ADDED FOR ALL OTHER CC MODELS.
C
C WHICH IS USED TO COMPUTE THE RLE ENERGY APPROXIMANT.

      SUBROUTINE FRMTAU_EXTRAP(Z,ETOT,ET2,TAU,ERLE,IORDER,SINGLE,RLETYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SINGLE
      INTEGER RLETYP
      DIMENSION ETOT(IORDER),ET2(IORDER),TAU(IORDER)
#ifdef _ASSERT
      if (iorder.lt.0) then
         print *, '@FRMTAU_EXTRAP: Assertion failed.'
         print *, '                iorder = ',iorder
         call errex
      end if
#endif /* _ASSERT */
      IF (RLETYP.EQ.1) THEN
         ERLE = 0.0
         DO I = 1, IORDER
            ERLE = ERLE + TAU(I)*ET2(IORDER+1-I)
         END DO
      ELSE
         ERLE = Z
         DO I = 1, IORDER
            ERLE = ERLE + TAU(I)*(ET2(IORDER+1-I)-Z)
         END DO
      END IF
      RETURN
      END 
