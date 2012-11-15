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
      subroutine apply_den4(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Each element of the array source(i,j,k,l) is divided by the 
c   corresponding element of the array target(i,j,k,l).  
c   Format of command is :
c
c   execute apply_den4 source target
c
c   There is no resriction on the type of array of source or target. 
c
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      if (nindex .ne. nindex2) then
         print *,'Error: Both arrays in apply_den4 must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex,' indices.'
         print *,'Second array has ',nindex2,' indices.'
         call abort_job()
      endif

      if (nindex .ne. 4) then 
         print *,'Arrays in apply_den4 must have 4 indices' 
         call abort_job()
      endif

      call do_div4(x,x2,
     *                      bval(1),eval(1),bval(2),eval(2),
     *                      bval(3),eval(3),bval(4),eval(4), 
     *                      bval2(1),eval2(1),bval2(2),eval2(2),
     *                      bval2(3),eval2(3),bval2(4),eval2(4)) 

      return
      end

      subroutine do_div4(x,y,
     *                   a1,a2,b1,b2,c1,c2,d1,d2,
     *                   i1,i2,j1,j2,k1,k2,l1,l2)
      implicit none 
      include 'interpreter.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer i1,i2,j1,j2,k1,k2,l1,l2
      double precision x(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision y(i1:i2,j1:j2,k1:k2,l1:l2)

      integer i,j,k,l 
      double precision val, div

      do i = i1, i2 
      do j = j1, j2 
      do k = k1, k2 
      do l = l1, l2 
         div = y(i,j,k,l) 
         val = x(i,j,k,l)/div   
         x(i,j,k,l) = val 
      enddo 
      enddo 
      enddo 
      enddo 

      return 
      end 

