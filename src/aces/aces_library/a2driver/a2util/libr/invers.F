
C THIS ARRAY CALCUATES THE INVERSE OF A GIVEN VECTOR A OF LENGTH N
C    A(I) = ONE/A(I)

      SUBROUTINE INVERS(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N)
      if (n.lt.1) return
      DO I = 1, N
         A(I) = 1.0d0 / A(I)
      END DO
      RETURN
      END
