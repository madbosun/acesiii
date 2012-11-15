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
      subroutine place_one_dea_7(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute place_one_dea_7 array scalar
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
         print *,'Error: place_one_dea_7 at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: place_one_dea_7 at line ',current_line,' not ',
     *           'called with scalar in 2nd arg.'
         call abort_job()      
      endif

      call doplaceonedea7(x, x2,
     *             bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine doplaceonedea7(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
c
      integer a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval   

      write(6,*) ' RANGE :', a1, a2, b1, b2, sval  
      write(6,*) ' alpha_occupied: ', nalpha_occupied  
      write(6,*) ' beta_occupied: ', nbeta_occupied  
      do a = a1, a2 
      do b = b1, b2 
        if(b .eq. nbeta_occupied+1 
     *        .and. a .eq. nalpha_occupied+3)then
           array1(a,b) = sval  
           write(6,*) 'C VECTOR ', a, b, array1(a,b)  
        endif
      enddo  
      enddo  
c
      return
      end

