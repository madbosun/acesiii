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
      subroutine return_h1(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c
c   return_h1 must be executed with the following syntax:
c   execute return_h1 h1(mu,nu)
c
c--------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      integer nsend

c---------------------------------------------------------------------------
c   Compute data block of integrals.
c---------------------------------------------------------------------------

      if (intpkg .eq.flocke_package) then
         call comp_return_h1(x, bval(1), eval(1), 
     *                       bval(2), eval(2),nsend)
      else
        write(6,*) ' MUST USE ERD IN RETURN_H1!!! ' 
      endif 

      return
      end
