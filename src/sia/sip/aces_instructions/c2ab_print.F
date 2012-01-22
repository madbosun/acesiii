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
      subroutine c2ab_print(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute c2ab_print array scalar
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      if (nindex .ne. 2) then
         print *,'Error: c2ab_print at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: c2ab_print at line ',current_line,' not ',
     *           'called with scalar in 2nd arg.'
         call abort_job()      
      endif

      call doc2ab_print(x, x2,
     *             bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine doc2ab_print(x,
     *                      i1,i2,j1,j2,k1,k2,l1,l2)
      implicit none 
      include 'interpreter.h'
      include 'parallel_info.h'

      integer i1,i2,j1,j2,k1,k2,l1,l2
      double precision x(i1:i2,j1:j2,k1:k2,l1:l2)

      integer a,b,c,d 
      double precision stemp
c
c      if (me .eq. 0)then
c
      do b = j1, j2
      do d = l1, l2
      do a = i1, i2
      do c = k1, k2
c
      stemp=dabs(x(a,b,c,d))
c
      if(stemp .gt. 0.05d-1)then
      write(6,*) ' C2AB VECTOR FOR I,J = ',b,d,',A,B = ',a,c,' = ',
     *x(a,b,c,d)
      endif
c
      enddo
      enddo
      enddo
      enddo
c
c
c      endif
c
c
      return 
      end 

