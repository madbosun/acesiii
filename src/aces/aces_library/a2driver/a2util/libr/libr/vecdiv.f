
C COMPUTES C(I)=A(I)/B(I) FOR FIRST N ELEMENTS OF VECTORS A AND B.

      SUBROUTINE VECDIV(A,B,C,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N),C(N)
      if (n.lt.1) return
      DO I = 1, N
         C(I) = A(I) / B(I)
      END DO
      RETURN
      END
