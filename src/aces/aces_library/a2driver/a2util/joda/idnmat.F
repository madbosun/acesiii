      SUBROUTINE IDNMAT(A,N,NTIME)
C
C SUBROUTINE RETURNS THE N X N IDENTITY MATRIX.  NTIME UTILITY IS
C  USED IN SYMUNQ.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision a(*)
      if (n.eq.3) then
         a(1) = 1.d0
         a(2) = 0.d0
         a(3) = 0.d0
         a(4) = 0.d0
         a(5) = 1.d0
         a(6) = 0.d0
         a(7) = 0.d0
         a(8) = 0.d0
         a(9) = 1.d0
      else
         do i = 1, n*n
            a(i) = 0.d0
         end do
         ndx = 1
         do i = 1, n
            a(ndx) = 1.d0
            ndx = ndx + n + 1
         end do
      end if
      NTIME=NTIME+1
      RETURN
      END
