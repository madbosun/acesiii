
C THIS ROUTINE COMPUTES AND WRITES OUT THE NORM OF a, WHICH IS AN
C VECTOR OF LENGTH n.

      subroutine znorm(a,n)
      implicit double precision (a-h,o-z)
      dimension a(n)
      z = 0.0d0
      if (n.gt.0) then
         do i = 1, n
            z = z + ( a(i) * a(i) )
         end do
      end if
      write(*,'(F20.10)') z
      return
      end
