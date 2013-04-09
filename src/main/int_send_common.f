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
      include 'scflindep.h'

      integer ierr, len
      integer*8 ixx, c_loc64

      common /flags/ iflags
      integer iflags(100) 

      integer ECP, ntot
      double precision ecp_erd
      common /FECP/ntot, ecp_erd(2500,2500)

c---------------------------------------------------------------------------
c   Broadcast int_gen_parms common block.
c---------------------------------------------------------------------------

      ixx = 1 
      len = (c_loc64(last, ixx, intsize) - 
     *       c_loc64(memptr, ixx, intsize)) / intsize
      call mpi_bcast(memptr, len, mpi_integer, 0, 
     *               mpi_comm_world, ierr)

c---------------------------------------------------------------------------
c  Broadcast fmo labeled common block.
c---------------------------------------------------------------------------

      call mpi_bcast(nfmo, 1001, mpi_integer, 0, 
     *               mpi_comm_world, ierr) 

c---------------------------------------------------------------------------
c  Broadcast iflags labeled common block.
c---------------------------------------------------------------------------

      call mpi_bcast(iflags, 101, mpi_integer, 0, 
     *               mpi_comm_world, ierr) 

c---------------------------------------------------------------------------
c  Broadcast ecp labeled common block.
c---------------------------------------------------------------------------

      call mpi_bcast(ntot, 1, mpi_integer, 0, 
     *               mpi_comm_world, ierr) 

      len = 2500*2500+1 
      call mpi_bcast(ecp_erd, len, mpi_double_precision, 0, 
     *               mpi_comm_world, ierr) 

c---------------------------------------------------------------------------
c     basis linear dependency threshold for SCF calculation
c---------------------------------------------------------------------------

      call mpi_bcast(lindep_tol,1,mpi_integer,0,mpi_comm_world, ierr)

c---------------------------------------------------------------------------
c     Integral threshold for SCF calculation
c---------------------------------------------------------------------------

      call mpi_bcast(intgrl_tol,1,mpi_integer,0,mpi_comm_world, ierr)
c---------------------------------------------------------------------------
c---------------------------------------------------------------------------
 
 
      return
      end
