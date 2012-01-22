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
      subroutine removeoo_dd(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      if (dryrun) return

       call retnooodd(x,x2,bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine retnooodd(array1,array2,a1,a2,b1,b2)

c---------------------------------------------------------------------------
c
c
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'

      integer i, a, b, c, d, a1, a2, b1, b2, c1, c2, d1, d2
      double precision array1(a1:a2,b1:b2)
      double precision array2(a1:a2,b1:b2)
      integer nzero, zero(100), index1, index2, index3, index4
      integer nszero, ndzero, szero(100), dzero(100)
      common /SINDEX/index1, index2, index3, index4
c
c     Determine range of doubly occupied indices to be removed.
c     ---------------------------------------------------------

      b      = 0
      ndzero = 0
      if (nalpha_occupied .gt. nbeta_occupied) then
         do a = 1, nbeta_occupied
            ndzero = ndzero + 1
            b = b + 1
            dzero(b) = a
         enddo
      endif
      if (nbeta_occupied .gt. nalpha_occupied) then
         do a = 1, nalpha_occupied
            ndzero = ndzero + 1
            b = b + 1
            dzero(b) = a
         enddo
      endif
c
c     Initialize output array to input array.
c     ---------------------------------------

      do a = a1, a2
      do b = b1, b2

         array2(a,b) = array1(a,b)

      enddo ! b
      enddo ! a
c
c     Reduce first index.
c     -------------------

      do i = 1,  ndzero
      do a = a1, a2
         if (dzero(i) .eq. a) then
         do b = b1, b2

            array2(a,b) = 0.0

         enddo ! b
         endif
      enddo ! a
      enddo ! i

100   continue
c
      return
      end
