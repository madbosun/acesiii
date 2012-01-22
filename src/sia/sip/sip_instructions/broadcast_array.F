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
      subroutine broadcast_array(x, nindex, type, bval,
     *                              eval, bdim, edim)
c---------------------------------------------------------------------------
c   Broadcasts the static array in the "result" field of the instruction
c   to each processor in the company.  The copy of the array on proc 0
c   is the data that is broadcast.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'parallel_info.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      integer i, n, company_comm, ierr
      integer pst_get_company_comm

c----------------------------------------------------------------------------
c   Broadcast the data.
c----------------------------------------------------------------------------

      n = 1
      do i = 1, nindex
         n = n * (edim(i)-bdim(i)+1)
      enddo

      company_comm = pst_get_company_comm(me)
      call mpi_bcast(x, n, MPI_DOUBLE_PRECISION, 0,
     *               company_comm, ierr)
      return
      end
