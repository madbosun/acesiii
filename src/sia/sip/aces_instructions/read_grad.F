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
      subroutine read_grad(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c   Divides each sample of the block given by the array argument by the 
c       eps = f_alpha(i)+f_beta(i)-f_alpha(a)-f_beta(a)
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      if (nindex .ne. 2) then
         print *,'Error in read_grad: ',
     *      'Gradient must be a 2-index array!.'
         print *,'array has ',nindex,' indices.'
         call abort_job()
      endif

      call add_grad(x,bval(1),eval(1),bval(2),eval(2))
      
      return
      end


      subroutine add_grad(x,iatoms,iatome,icomponents,icomponente)
      implicit none
      include 'parallel_info.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'

      integer i,j,igrad  
      integer iatoms,iatome,icomponents,icomponente  
      double precision x(iatoms:iatome,icomponents:icomponente)

      if (me .ne. 0) return  

      igrad = (iatoms-1)*3 + icomponents 

      gradient_data(igrad) = gradient_data(igrad)+ x(iatoms,icomponents) 

      write(6,*) ' GRAD :', iatoms,icomponents,igrad 

33    format(3F16.8)

      return
      end


