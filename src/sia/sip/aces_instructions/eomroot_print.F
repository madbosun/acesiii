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
      subroutine eomroot_print(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute eomroot_print array scalar
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
         print *,'Error: eomroot_print at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: eomroot_print at line ',current_line,' not ',
     *           'called with scalar in 2nd arg.'
         call abort_job()      
      endif

      call doeomroot_print(x, x2,
     *             bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine doeomroot_print(array1,sval,a1,a2,b1,b2) 

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

      if (me .eq. 0)then
      if (sval .eq. 1) then 
      write(6,*) ' ' 
      write(6,*) ' ---------------------------------------- ' 
      write(6,*) '             EOM-CCSD energies ' 
      write(6,*) '  root       Hartree        eV '  
      write(6,*) ' ---------------------------------------- ' 
      endif 
c
      do a = a1, a2 
         stemp=array1(a,a)
         write(6,100) a, stemp, stemp*27.2113957 
      enddo  

      if (sval .eq. 12) then 
      write(6,*) ' ' 
      write(6,*) ' ---------------------------------------- ' 
      write(6,*) ' ---------------------------------------- ' 
      write(6,*) ' ' 
      endif 
c
      endif

100   format(I4,5x,F12.8,5x,F12.8)  
c
      return
      end

