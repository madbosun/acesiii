
C THIS ROUTINE ACCEPTS A FLOATING POINT VECTOR (VEC) AND
C RETURNS VEC(I)=1.0/SQRT(VEC(I)) FOR THE FIRST n
C ELEMENTS FOR STRIDE OFFSET inc.
c BEWARE - Do not use this (yet) if inc is negative!

      SUBROUTINE INVSQT(VEC,inc,n)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION VEC(*)
c old
c      IOFF = 1
c      DO I = 1, ICOUNT
c         VEC(IOFF) = 1.0d0 / SQRT(VEC(IOFF))
c         IOFF = IOFF + ISTRID
c      END DO
c new
      if (n.lt.1) return
      if (inc.eq.1) then
         do i = 1, n
            vec(i) = 1.0d0 / sqrt(vec(i))
         end do
      else
         itmp = 1 + inc*(n-1)
         do i = 1, itmp, inc
            vec(i) = 1.0d0 / sqrt(vec(i))
         end do
      end if
c end
      RETURN
      END
