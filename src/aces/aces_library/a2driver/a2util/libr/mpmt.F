
C SUBROUTINE ADDS SQUARE MATRIX X AND ITS TRANSPOSE AND PUTS RESULT IN
C ORIGINAL MATRIX.

      SUBROUTINE MPMT(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,N)
c old
c      DO I = 1, N
c         DO J = 1, I
c            Z      = X(I,J) + X(J,I)
c            X(I,J) = Z
c            X(J,I) = Z
c         END DO
c      END DO
c new
      if (n.lt.1) return
      if (n.eq.1) then
         x(1,1) = 2.0d0 * x(1,1)
      else
         do j = 1, n
            do i = j, n
               x(i,j) = x(i,j) + x(j,i)
            end do
         end do
         do j = 1, n
            do i = 1, j
               x(i,j) = x(j,i)
            end do
         end do
      end if
c end
      RETURN
      END
