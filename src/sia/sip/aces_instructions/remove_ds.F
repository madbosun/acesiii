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
      subroutine remove_ds(x,nindex, type, bval,
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

       call retnods(x,x2,bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine retnods(array1,array2,a1,a2,b1,b2)

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

c     Initialize output array to input array.
c     ---------------------------------------

      do a = a1, a2
      do b = b1, b2

         array2(a,b) = array1(a,b)

      enddo ! b
      enddo ! a
c
c     Determine range of singly occupied indices to be removed.
c     ---------------------------------------------------------

      b      = 0
      nszero = 0
      if (nalpha_occupied .gt. nbeta_occupied) then
         nszero = nalpha_occupied - nbeta_occupied
         do a = nalpha_occupied, nbeta_occupied, -1
            b = b + 1
            szero(b) = a
         enddo
      endif
      if (nbeta_occupied .gt. nalpha_occupied) then
         nszero = nbeta_occupied - nalpha_occupied
         do a = nbeta_occupied, nalpha_occupied, -1
            b = b + 1
            szero(b) = a
         enddo
      endif

      if (nalpha_occupied .eq. nbeta_occupied) go to 100
c
c     Reduce second index.
c     -------------------

      do i = 1,  nszero
      do b = b1, b2
         if (szero(i) .eq. b) then
         do a = a1, a2

            array2(a,b) = 0.0

         enddo ! b
         endif
      enddo ! a
      enddo ! i
c
c     Determine range of doubly occupied indices to be removed.
c     ---------------------------------------------------------

      b      = 0
      ndzero = 0
      if (nalpha_occupied .gt. nbeta_occupied) then
         ndzero = nbeta_occupied
         do a = 1, nbeta_occupied
            b = b + 1
            dzero(b) = a
         enddo
      endif
      if (nbeta_occupied .gt. nalpha_occupied) then
         ndzero = nalpha_occupied
         do a = 1, nalpha_occupied
            b = b + 1
            dzero(b) = a
         enddo
      endif
c
c     Reduce first index.
c     -------------------

      do i = 1,  ndzero
      do a = a1, a2
         if (dzero(i) .eq. a) then
         do b = b1, b2

            array2(a,b) = 0.0

         enddo ! a
         endif
      enddo ! b
      enddo ! i

100   continue
c
      return
      end

