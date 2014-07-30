C
      SUBROUTINE SET_2UNIT_MATRIX(A, LEN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LEN, LEN)
C
      if (len.lt.1) return
 
      CALL ZERO(A, LEN*LEN)
C
      DO I = 1, LEN
         A(I, I) = 1.0D0
      END DO
C
      RETURN
      END
