C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine symm_force_a(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c The diagonal elements of a four index array are extracted. 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      if (nindex .eq. 4) then
         call dosymmforce(x, bval(1),eval(1),type(1),
     *                   bval(2),eval(2),type(2),
     *                   bval(3),eval(3),type(3),
     *                   bval(4),eval(4),type(4), 
     *                   epsilon,epsilonb)
      else
         print *,'Error in symm_force_a: Result array ',
     *      'must have 4 indices.'
         print *,'array has ',nindex,' indices.'
         call abort_job()
      endif
      
      return
      end

      subroutine dosymmforce(x, a1,a2,atype,
     *                      i1,i2,itype,b1,b2,btype,
     *                      j1,j2, jtype, epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j
      double precision val
c
      if (i1 .eq. 1 .and. i2 .eq. 1) then 
         do a = a1,a2
         do b = b1,b2
            do i = i1,i2
            do j = j1,j2
               if (a .lt. b)then
                  val = x(a,i,b,j)
                  x(b,i,a,j)= val*(-1.0d0)
               endif
               if (a .eq. b)then
                  val = 0.0 
                  x(b,i,a,j)= val*(-1.0d0)
               endif
            enddo
            enddo
         enddo
         enddo
         return 
      endif 
c
      do b = b1,b2
      do a = a1,a2
         do j = j1,j2
         do i = i1,i2
            if (a .lt. b)then
               val = x(a,i,b,j)
               x(b,i,a,j)= val*(-1.0d0)
            endif
         enddo
         enddo
      enddo
      enddo
c
      val = 0.0d0 
      do b = b1,b2
      do a = a1,a2
         do j = j1,j2
         do i = i1,i2
         if(i .eq. j .or. a .eq. b)then
               x(a,i,b,j)= val
         endif
         enddo
         enddo
      enddo
      enddo
c
      return
      end

