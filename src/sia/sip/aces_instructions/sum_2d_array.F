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
      subroutine Sum_2D_Array(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                       x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute Sum_2D_Array array1 sum
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
         print *,'Error: return_sval at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: return_sval at line ',current_line,' not ',
     *           'called with scalar in 2nd arg. nindex2 =', nindex2 
         call abort_job()
      endif

      call return_sum(x, x2,
     *                bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine return_sum(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      integer a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval   

      sval = 0.0d+0 
      do a = a1, a2 
         do b = b1, b2 
           sval = Sval + array1(a,b) 

#ifdef _DEBUG_LVL0
           write(6,"(a,1x,2I4,1x,F20.10) ' SVAL ', a, b, sval 
#endif
      enddo ! b 
      enddo ! a 
      sval = array1(a2,b2)
c
      return
      end

