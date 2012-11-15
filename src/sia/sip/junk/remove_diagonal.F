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
      subroutine remove_diagonal(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c
c   The diagonal elements of an array(Array1) are removed and the array 
c   without the diagonal elements is put into another array(Array2). 
c
c   The instruction is executed as follows:
c   
c   execute remove_diagonal Array1 Array2
c
c   where Array1 is the full input matrix, usually the fock matrix.  
c         Array2 is the matrix with zero diagonal elements.  
c
c   Both Array1 and Array2 must be "static" arrays.
c--------------------------------------------------------------------------

      implicit none
      include 'trace.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      integer n

      if (dryrun) return

      n = edim(1) - bdim(1) + 1 
      if (nindex .eq. 2) call rdiag2(x,x2,n) 

      return 
      end 

      subroutine rdiag2(array1,array2,n) 

c---------------------------------------------------------------------------
c
c The diagonal elements of the array1 are removed and the output put into 
c array2. The diagonal elements of the Fock matrix are also put into
c the epsilon array.  
c 
c---------------------------------------------------------------------------

      implicit none
      include 'epsilon.h'
      integer n 
      double precision array1(n,n)
      double precision array2(n,n)
      integer a, b

      do b = 1, n
      do a = 1, n
         if (a .ne. b) array2(a,b) = array1(a,b)
         if (a .eq. b) array2(a,b) = 0.0  
c        write(6,*) a, b, array1(a,b), array2(a,b) 
      enddo
      enddo

      return
      end

