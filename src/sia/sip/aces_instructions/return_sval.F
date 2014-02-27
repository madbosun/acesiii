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
      subroutine return_sval(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute return_sval array scalar
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

      if (nindex .gt. 2) then
         print *,'Error: return_sval at line ',current_line
         print *,'First arg. must be a block of 1 or 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: return_sval at line ',current_line,' not ',
     *           'called with scalar in 2nd arg. nindex2 =', nindex2 
         call abort_job()
      endif

      if (nindex .eq. 2) 
     *    call doreturn2(x, x2,
     *                  bval(1),eval(1),bval(2),eval(2))

      if (nindex .eq. 1) 
     *    call doreturn1(x, x2,
     *                  bval(1),eval(1))

      return 
      end 

      subroutine doreturn2(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      integer a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval   

      sval = 0.0d+0 
c     write(6,*) ' RANGE :', a1, a2, b1, b2 
c      do a = a1, a2 
c      do b = b1, b2 
c           sval = array1(a,b) 
c          write(6,*) ' SVAL ', a, b, sval 
c      enddo ! b 
c      enddo ! a 
c     write(6,*) ' Maximum density difference :', val 
      sval = array1(a2,b2)
c
      return
      end

      subroutine doreturn1(array1,sval,a1,a2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      integer a, a1, a2
      double precision array1(a1:a2), sval   

      sval = 0.0d+0 
c     write(6,*) ' RANGE :', a1, a2, b1, b2 
c      do a = a1, a2 
c      do b = b1, b2 
c           sval = array1(a,b) 
c          write(6,*) ' SVAL ', a, b, sval 
c      enddo ! b 
c      enddo ! a 
c     write(6,*) ' Maximum density difference :', val 
      sval = array1(a2)
c
      return
      end

