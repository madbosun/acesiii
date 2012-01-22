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
      subroutine c1b_print(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute c1b_print array scalar
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
         print *,'Error: c1b_print at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: c1b_print at line ',current_line,' not ',
     *           'called with scalar in 2nd arg.'
         call abort_job()      
      endif

      call doc1b_print(x, x2,
     *             bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine doc1b_print(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   

c      if (me .eq. 0)then
c
      do b = b1, b2 
      do a = a1, a2 
      stemp=dabs(array1(a,b))
c
      if(stemp .gt. 0.05d-1)then
      write(6,*) ' C1B VECTOR FOR I = ',b,',A = ',a,' = ',array1(a,b)
      endif

C      if(stemp .gt. 0.0001d-1)then 
C      write(6,999) ' C1B VECTOR FOR I = ',b,',A = ',a,' = ',array1(a,b)  
C      endif
C 999  FORMAT (A20,I5,A5,I5,A3,F14.8)

c
      enddo  
      enddo  
c
c
c      endif
c
      return
      end

