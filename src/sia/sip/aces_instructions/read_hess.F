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
      subroutine read_hess(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c     The 4-index array (hessian at present) is read in and summed 
c     into an array in a common block. 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'int_gen_parms.h'
      include 'hess.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      double precision xx(1)
#ifdef ALTIX
      pointer (dptr, xx)
#else
      common xx
#endif

      if (dryrun) return
      write_hess = .true.

#ifdef ALTIX
      dptr = dshptr
#endif

      if (nindex .ne. 4) then
         print *,'Error in read_hess: ',
     *      'Hessian must be a 4-index array!.'
         print *,'array has ',nindex,' indices at line ',
     *            current_line,'.'
         call abort_job()
      endif

      call add_hess(x,bval(1),eval(1),bval(2),eval(2),
     *                        bval(3),eval(3),bval(4),eval(4),
     *                        xx(ix_hess))
      
      return
      end


      subroutine add_hess(x,iatoms,iatome,icomponents,icomponente,
     *                      jatoms,jatome,jcomponents,jcomponente, 
     *                      hess) 
      implicit none
      include 'parallel_info.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'

      integer i,j,ihess,jhess 
      integer iatoms,iatome,icomponents,icomponente  
      integer jatoms,jatome,jcomponents,jcomponente  
      double precision x(iatoms:iatome,icomponents:icomponente,
     *                   jatoms:jatome,jcomponents:jcomponente)
      double precision hess(3*ncenters,3*ncenters) 


      if (dbg) then
         write(6,*) ' HCOORD :', iatoms,icomponents,jatoms,
     *                    jcomponents
         write(6,*) ' Hessian from Loops ',x(iatoms,icomponents, 
     *                                    jatoms,jcomponents)  
      endif

      if (me .ne. 0) return  

      ihess = (iatoms-1)*3 + icomponents 
      jhess = (jatoms-1)*3 + jcomponents 

      hess(jhess,ihess) = hess(jhess,ihess) 
     *                  + x(iatoms,icomponents,jatoms,jcomponents) 

c     if (ihess .ne. jhess) 
c      hess(jhess,ihess) = hess(jhess,ihess) 
c     *                  + x(iatoms,icomponents,jatoms,jcomponents) 

33    format(3F16.8)

      return
      end


