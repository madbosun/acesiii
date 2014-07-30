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
      subroutine read_density(watom,
     *                 scr, maxblk, iscr, coords, coeffs, alphas,
     *                 ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 nai, kin, ovl,
     *                 fa,fb)
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'parallel_info.h'

      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer a,b,c,d
      integer iatom, n_basis
      double precision watom

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst
      integer ierr
      integer nc1, nc2, nd1, nd2

      integer imin, zmin, iblk, zblk

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision x3,y3,z3
      double precision x4,y4,z4

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision nai(nc1:nc2,nd1:nd2)
      double precision kin(nc1:nc2,nd1:nd2)
      double precision ovl(nc1:nc2,nd1:nd2)
      double precision H0T(nc1:nc2,nd1:nd2)

      double precision ca(nc1:nc2,nc1:nc2)
      double precision cb(nc1:nc2,nc1:nc2)

      double precision fa(nc1:nc2,nc1:nc2)
      double precision fb(nc1:nc2,nc1:nc2)
      integer map(nc1:nc2)
      integer umap(nc1:nc2)
      integer beg_anfps(max_shells)
      integer end_anfps(max_shells)
      double precision scr(*)
      integer iscr(*)

      integer ccbeg(*), ccend(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      integer*8 ccbeg_pack64(max_dim_coeff), ccend_pack64(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff),
     *                 pcoeff_pack(max_dim_coeff)
      integer*8 arg64(25)
      integer n_max, natoms
      logical aExist, bExist, breturn
      integer Ierror

      common /Imax_com/sz_max(max_shells,max_shells), delta
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter

      common /orbread/breturn

      save alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

      open(66,file='summary.out')

c     call mpi_comm_rank(mpi_comm_world, me, ierr)

      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

c-----------------------------------------------------------------------
c   Check if a guess file exists. If yes READ the orbitals, 
c   construct the density and finish.  
c-----------------------------------------------------------------------

      if (breturn) return

      nbasis = nalpha_occupied + nalpha_virtual
101   Format(5x,'Something wring with read orbitals. N read N correct',
     *       2I8)

      INQUIRE(file='ca.data',exist=aExist)
      INQUIRE(file='cb.data',exist=bExist)

      if (aexist .and. me .ne. 0) return
      if (bexist .and. me .ne. 0) return

      if (aexist) then
      if (bexist) then
      else

         if (me .eq. 0 .and. .not. breturn) then
            OPEN (135, File = 'ca.data',
     *            Status = 'Old', Iostat = Ierror)

            rewind 135
            read(135,*) n_basis
            write(66,*) ' Reading alpha orbitals'
                        if (n_basis .ne. nbasis)
     *                  write(66,101) n_basis, nbasis
            do a = 1, n_basis
            do b = 1, n_basis
               read(135,*) ca(a,b)
               cb(a,b) = ca(a,b)
            enddo
            enddo

            close (135)

            write(66,*) ' Computing the RHF HF density'
            call hfdensity(ca,cb,FA,FB,nbasis,nalpha_occupied,
     *                     nbeta_occupied)
            breturn = .true.
            watom   = 1.0
            return

         endif

      endif ! aexist  
      endif ! bexist 

      if (aexist .and. bexist) then

         if (me .eq. 0 .and. .not. breturn) then
            OPEN (135, File = 'ca.data',
     *            Status = 'Old', Iostat = Ierror)

            rewind 135
            read(135,*) n_basis
            write(66,*) ' Reading alpha orbitals'
                        if (n_basis .ne. nbasis)
     *                  write(66,101) n_basis, nbasis
            do a = 1, n_basis
            do b = 1, n_basis
               read(135,*) ca(a,b)
            enddo
            enddo

            close (135)
            OPEN (357, File = 'cb.data',
     *            Status = 'Old', Iostat = Ierror)

            rewind 357
            read(357,*) n_basis
            write(66,*) ' Reading beta orbitals'
                        if (n_basis .ne. nbasis)
     *                  write(66,101) n_basis, nbasis
            do a = 1, n_basis
            do b = 1, n_basis
               read(357,*) cb(a,b)
            enddo
            enddo
            close(357)

            write(66,*) ' Computing the UHF density'
            call hfdensity(ca,cb,FA,FB,nbasis,nalpha_occupied,
     *                     nbeta_occupied)
            breturn = .true.
            watom   = 1.0
            return

         endif

      endif ! aexist and bexist 

      if (breturn) return
c-----------------------------------------------------------------------
c   Done initial guess from file  
c-----------------------------------------------------------------------

      end

