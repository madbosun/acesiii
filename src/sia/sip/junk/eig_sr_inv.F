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
      subroutine eig_sr_inv(x, nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Calculates the eigenvalues and eigenvectors of a 2-d square matrix.
c   The matrix is also diagonalized on output and inverted. 
c
c   The eigenvalues are stored in a common block.  
c   The instruction is executed as follows:
c   
c   execute eig Array1 Array2
c
c   where Array1 is the input matrix that is to be diagonalized.
c         Array2 is the matrix whose columns are the eigenvectors of Array1.
c
c   Both Array1 and Array2 must be "static" arrays.
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

      integer n, ierr, junk

      ierr = 0
      if (nindex .ne. 2) then
         ierr = ierr + 1
         print *,'Error in eig_sr_inv, line ',current_line,
     *           ' Number of array',
     *           ' indices must be 2, actual number is ',nindex
      endif

      if (nindex2 .ne. 2) then
         ierr = ierr + 1
         print *,'Error in eig_sr_inv, line ',current_line,
     *           ' Number of array',
     *           ' indices must be 2, actual number is ',nindex2
      endif

      if (ierr .ne. 0) call abort_job()

      n = edim(1) - bdim(1) + 1

      call diag(x, x2, junk, n, 1, 1, 1)

      return
      end
