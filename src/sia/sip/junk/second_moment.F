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
      subroutine second_moment(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute second_moment array scalar
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      integer itype, nsend

      if (nindex .ne. 2) then
         print *,'Error: second_moment at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: second_moment at line ',current_line,' not ',
     *           'called with scalar in 2nd arg.'
         call abort_job()      
      endif

      itype = 0
      if ( (abs (x2 - 1.0)) .lt. 0.001 ) itype = 1
      if ( (abs (x2 - 2.0)) .lt. 0.001 ) itype = 2
      if ( (abs (x2 - 3.0)) .lt. 0.001 ) itype = 3
      if ( (abs (x2 - 4.0)) .lt. 0.001 ) itype = 4
      if ( (abs (x2 - 5.0)) .lt. 0.001 ) itype = 5
      if ( (abs (x2 - 6.0)) .lt. 0.001 ) itype = 6


c---------------------------------------------------------------------------
c   Compute data block of integrals.
c---------------------------------------------------------------------------

      if (intpkg .eq.flocke_package) then
         call comp_second_moment(itype,x, bval(1), eval(1), 
     *                           bval(2), eval(2), nsend)
      else
        write(6,*) ' MUST USE ERD IN SECOND_MOMENT!!! ' 
        call abort_job()
      endif 

      return
      end
