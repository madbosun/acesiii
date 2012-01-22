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
      subroutine return_diagonal4(x, nindex, type, bval,
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
         call ret_diag4(x, bval(1),eval(1),type(1),
     *                   bval(2),eval(2),type(2),
     *                   bval(3),eval(3),type(3),
     *                   bval(4),eval(4),type(4), 
     *                   epsilon,epsilonb)
      elseif (nindex .eq. 2) then
         call lret_diag2(x, bval(1),eval(1),type(1),
     *                   bval(2),eval(2),type(2),
     *                   epsilon,epsilonb)
      else
         print *,'Error in return_diagonal_4: Result array ',
     *      'must have 2 or 4 indices.'
         print *,'array has ',nindex,' indices.'
         call abort_job()
      endif
      
      return
      end

      subroutine ret_diag4(x, a1,a2,atype,
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

      do j = j1,j2
      do b = b1,b2
         do i = i1, i2
         do a = a1,a2
            val = 0.0d0 
            if ((atype .eq. itype) .and. (i .eq. a) .and. 
     *          (btype .eq. jtype) .and. (j .eq. b)) then  
               val =  x(a,i,b,j)
            else 
               val = 0.0d0 
            endif 
            x(a,i,b,j) = val 
         enddo
         enddo
      enddo
      enddo

      return
      end

      subroutine lret_diag2(x, a1,a2,atype,i1,i2,itype,
     *                      epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j
      double precision val

      do i = i1, i2
      do a = a1 ,a2
         val = 0.0d0 
         if ((atype .eq. itype) .and. (i .eq. a)) then 
            val =  x(a,i)
         else 
            val = 0.0d0 
         endif 
         x(a,i) = val 
      enddo
      enddo

      return
      end

