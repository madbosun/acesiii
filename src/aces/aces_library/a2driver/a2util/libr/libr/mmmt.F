
C SUBROUTINE FORMS X-X(transpose) AND RETURNS RESULT IN X.

      SUBROUTINE MMMT(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,N)
c old
c      DO I = 1, N
c         DO J = 1, I
c            dtmp   = X(I,J) - X(J,I)
c            X(I,J) =  dtmp
c            X(J,I) = -dtmp
c         END DO
c      END DO
c new
      if (n.lt.1) return
      if (n.eq.1) then
         x(1,1) = 0.0d0
      else
         do j = 1, n
            do i = j, n
               x(i,j) = x(i,j) - x(j,i)
            end do
         end do
         do j = 1, n
            do i = 1, j
               x(i,j) = -x(j,i)
            end do
         end do
      end if
c end
      RETURN
      END
