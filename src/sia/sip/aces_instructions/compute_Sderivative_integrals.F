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
      subroutine compute_Sderivative_integrals(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c   CHANGE COMMENT 
c   Multiplies the each sample of the block given by the array argument by 
c   diagonal eigenvalue matrix.  The eigenvalues matrix is not stored in
c   a 2-dimensional form, since it can be represented by a vector of
c   the diagonal elements.
c--------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'
      include 'hess.h'
      include 'trace.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (iptr, ibuf)
#else
      equivalence (ibuf, dbuf)
      common dbuf
#endif

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      integer der_flags(12)
      integer k

      if (dryrun) return
 
#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c-----------------------------------------------------------------------------
c   Calculate the integral blocks.
c-----------------------------------------------------------------------------

         do k = 1, 12
            der_flags(k) = 0
         enddo

c-----------------------------------------------------------------------------
c   Computation of Hessian contribution.
c-----------------------------------------------------------------------------

         call contract_erd_2der(bval(1),eval(1),bval(2),eval(2),
     *               bval(3),eval(3),bval(4),eval(4),
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               x, dbuf(ix_hess))

      go to 999
c-----------------------------------------------------------------------------
c   Test of second-derivatives.
c-----------------------------------------------------------------------------

         write(6,*) ' Computing Second-derivatives '
         call tester_MRT2(bval(1),eval(1),bval(2),eval(2),
     *               bval(3),eval(3),bval(4),bval(4),
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               x, dbuf(ix_hess))

c-----------------------------------------------------------------------------
c   Test of third-derivatives.
c-----------------------------------------------------------------------------

         write(6,*) ' Computing Third-derivatives '
         call tester_MRT3(bval(1),eval(1),bval(2),eval(2),
     *               bval(3),eval(3),bval(4),bval(4),
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               x)
c-----------------------------------------------------------------------------
c   Test of fourth-derivatives.
c-----------------------------------------------------------------------------

         write(6,*) ' Computing Fourth-derivatives '
         call tester_MRT4(bval(1),eval(1),bval(2),eval(2),
     *               bval(3),eval(3),bval(4),bval(4),
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               x)

999   continue

      return
      end

