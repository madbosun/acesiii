
c finds dsqrt(a*a+b*b) without overflow or destructive underflow

      double precision function pythag(a,b)
      implicit none

      double precision a,b
      double precision p,r,s,t,u,dtmp
      integer i

      p = max(abs(a),abs(b))

      if (p .eq. 0.0d0) then
         pythag = 0.0d0
         return
      end if

      r = min(abs(a),abs(b)) / p
      r = r*r

      i = 0
      do while (i.lt.100)
         t = 4.0d0 + r
         if (t .eq. 4.0d0) then
            pythag = p
            return
         end if
         s = r * (1/t)
         u = 1.0d0 + ( 2.0d0 * s )
         p = p * u
         dtmp = s * (1/u)
         r = r * dtmp * dtmp
         i = i + 1
      end do

c     At this point, the careful algorithm did not converge.
c     Screw it; just calculate the square root.

#ifdef _DEBUG
      write(*,*) '@PYTHAG: The careful algorithm did not converge.'
      write(*,*) '         Calculating by brute force:'
      write(*,*) '         pythag = sqrt( ',a,'**2 + ',b,'**2 )'
#endif
      pythag = sqrt( (a*a) + (b*b) )

      return
      end

