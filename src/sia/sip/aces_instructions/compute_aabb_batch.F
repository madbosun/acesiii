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
      subroutine compute_aabb_batch(x, nindex, type, bval,
     *                              eval, bdim, edim)
c-------------------------------------------------------------------------- 
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (iptr, ibuf)
#else
      common dbuf
      equivalence (ibuf(1), dbuf(1))
#endif

      integer*8 integral_scr
      integer maxblk , nsend

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Compute the block of integrals.
c---------------------------------------------------------------------------

         if (intpkg .eq. flocke_package) then
            integral_scr = dpkgscr
            maxblk  = zmax
         else if (intpkg .eq. gamess_package .or.
     *            intpkg .eq. gamess_derivative_package) then
            integral_scr = i0wrk
            maxblk = 15**4    ! max. angular momentum for GAMESS is 15.
         endif

      if (.not. dryrun) 
     *   call compute_aabb_integrals(bval(1),eval(1), bval(2),eval(2),
     *                       bval(3),eval(3),bval(4),eval(4),
     *                 dbuf(integral_scr),maxblk, ibuf(ipkgscr),
     *                 dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha), 
     *                 ibuf(iccbeg), ibuf(iccend), 
     *                 x, nsend)

      return
      end
