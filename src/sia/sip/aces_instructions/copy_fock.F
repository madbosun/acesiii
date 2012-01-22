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
      subroutine copy_fock(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Copies one array into a Fock array and copies the diagonal elements
c   into the corresponding eigenvalue array.
c
c   execute Array1 fock_a 
c         OR
c   execute Array1 fock_b
c 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'

      double precision x(*) 
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      integer i, n
 
      n = edim(1) - bdim(1) + 1 
      do i = 1, n*n
         x2(i) = x(i)
      enddo

c-------------------------------------------------------------------------
c   Copy diagonal elements into epsilon or epsilonb.
c-------------------------------------------------------------------------

      if (type(1) .eq. moaindex) then
         call copy_diag(x2, n, epsilon)
      else
         call copy_diag(x2, n, epsilonb)
      endif

      return
      end

      subroutine copy_diag(x, n, y)
      implicit none
      integer i, n
      double precision x(n,n), y(n)

      do i = 1, n
         y(i) = x(i,i)
      enddo
      return
      end
