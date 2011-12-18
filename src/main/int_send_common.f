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
      subroutine int_send_common()
c-------------------------------------------------------------------------
c   Broadcast the master's common block.
c------------------------------------------------------------------------- 
      implicit none

      include 'mpif.h'
      include 'proto_events.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'fmo.h'
      include 'shared_mem_data.h'

      integer ierr, len
      integer*8 ixx, c_loc64

c---------------------------------------------------------------------------
c   Broadcast int_gen_parms common block.
c---------------------------------------------------------------------------

      ixx = 1 
      len = (c_loc64(last, ixx, intsize) - 
     *       c_loc64(memptr, ixx, intsize)) / intsize
      call mpi_bcast(memptr, len, mpi_integer, 0, 
     *               mpi_comm_world, ierr)

c---------------------------------------------------------------------------
c   Broadcast ixalpha block, which is allocated using shared memory
c---------------------------------------------------------------------------
c      call mpi_bcast(shared_buf(shared_mem_offset),
c     *     max_shells, mpi_integer, 0, mpi_comm_world, ierr)

      call bcast_shared_mem()
c     Do we need a barrier here?
      call mpi_barrier(mpi_comm_world,ierr)
c      call protect_shared_mem()
c---------------------------------------------------------------------------
c  Broadcast fmo labeled common block.
c---------------------------------------------------------------------------

      call mpi_bcast(nfmo, 1001, mpi_integer, 0, 
     *               mpi_comm_world, ierr) 
 
      return
      end
