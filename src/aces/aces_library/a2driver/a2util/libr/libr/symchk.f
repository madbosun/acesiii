
C THIS ROUTINE CHECKS THE SYMMETRY OF AN NxN DOUBLE PRECISION 
C MATRIX (A).  FOR ALL ELEMENTS WHICH DIFFER BY MORE THAN 
C TOL WITH THEIR PARTNER [(I,J) AND (J,I)], THE INDICES AND
C VALUES ARE PRINTED OUT [I,J,A(I,J) AND A(J,I)].

      subroutine symchk(a,n,tol)
      implicit double precision (a-h,o-z)
      dimension a(n,n)
      if (n.lt.1) return
      do i = 1, n
         do j = 1, i-1
            if (dabs(a(i,j)-a(j,i)).gt.tol) then
               write(*,'(t3,2I5,2F10.7)') i,j,a(i,j),a(j,i)
            end if
         end do
      end do
      return
      end
