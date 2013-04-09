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
      subroutine aces_to_erd(nfps, iangular, nshells, ispherical,
     *                       erd_index, scalars) 
c-----------------------------------------------------------------------------
c   Using the shell angular momentum and basis function information, this
c   subroutine calculates an array which maps the ACES integral order to 
c   the ERD integral order.  This array may be used to re-order a block 
c   of integrals calculated by the ERD package into a format corresponding
c   to the ACES (VMOL-based) integrals.
c
c   Arguments:
c      nfps		Number of basis functions per shell.
c      iangular		Shell type (based on angular momentum) per shell.
c			(i. e. an S shell = 0, P = 1, etc.)
c      nshells		Number of shells contained in nfps and iangular.
c      ispherical       1 = spherical coordinates, 0 = Cartesian
c      erd_index	An array (output) of nbasis indices.  The ith
c			value of erd_index is the index in the ERD-based 
c        		system corresponding to index "i" in the ACES-based
c			system.
c      scalars          An array (output) of scale factors to adjust the
c                       ERD integrals to match VMOL integrals.
c-----------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'saved_data.h'
c     include 'int_gen_parms.h'

      integer n, nshells, ispherical
      integer nfps(nshells), iangular(nshells)
      integer erd_index(*)
      double precision scalars(*)

      integer istart, iend, ishell
      integer i, j, ierd
      integer nctr, nsh_coords
  
      integer smap(1)   
      integer pmap(3)
      integer dmap(5)
      integer fmap(7)
      integer gmap(9)
      integer hmap(11)
      integer imap(13)

      double precision d_scale(5)
      double precision f_scale(7)
      double precision g_scale(9)
      double precision h_scale(11)
      double precision i_scale(13)
      double precision xscale

      double precision x_scale(0:6,1:28) 

      common /flags/ iflags
      integer iflags(100)

      save smap, pmap, dmap, fmap, gmap, hmap, imap
      save d_scale, f_scale, g_scale, h_scale, i_scale

      data smap/1/
      data pmap/1,2,3/
      data dmap/4,3,1,5,2/
      data fmap/4,7,1,3,2,5,6/
      data gmap/4,8,6,3,1,9,2,5,7/
      data hmap/10,7,4,3,1,9,2,11,8,5,6/
      data imap/1,3,4,8,6,13,10,12,2,11,9,5,7/
    
      if (first_aces_to_erd) then
         first_aces_to_erd = .false.
         d_scale(1) = 2.d0 * dsqrt(3.d0)
         d_scale(2) = 1.d0
         d_scale(3) = 1.d0
         d_scale(4) = 2.d0
         d_scale(5) = 1.d0

         f_scale(1) = 2.d0 * dsqrt(10.d0)
         f_scale(2) = 2.d0 * dsqrt(10.d0)
         f_scale(3) = 2.d0 * dsqrt(15.d0)
         f_scale(4) = 2.d0 * dsqrt(6.d0)
         f_scale(5) = 1.d0
         f_scale(6) = 2.d0 * dsqrt(6.d0)
         f_scale(7) = 2.d0  

         g_scale(1) = 8.d0 * dsqrt(105.d0)
         g_scale(2) = 2.d0 * dsqrt(21.d0)
         g_scale(3) = 2.d0 * dsqrt(42.d0)
         g_scale(4) = 8.d0 * dsqrt(3.d0)
         g_scale(5) = 2.d0 * dsqrt(6.d0)
         g_scale(6) = 4.d0 * dsqrt(21.d0)
         g_scale(7) = 2.d0 * dsqrt(3.d0)
         g_scale(8) = 2.d0 * dsqrt(6.d0)
         g_scale(9) = 2.d0 * dsqrt(42.d0)

         h_scale(1) = 24.d0 * dsqrt(7.d0)
         h_scale(2) = 24.d0 * dsqrt(7.d0)
         h_scale(3) = 12.d0
         h_scale(4) = 24.d0 * dsqrt(6.d0)
         h_scale(5) = 2.d0 * dsqrt(3.d0)
         h_scale(6) = 8.d0 * dsqrt(30.d0)
         h_scale(7) = 8.d0 * dsqrt(3.d0)
         h_scale(8) = 24.d0 * dsqrt(6.d0)
         h_scale(9) = 24.d0 * dsqrt(105.d0)
         h_scale(10) = 8.d0 * dsqrt(30.d0)
         h_scale(11) = 6.d0

         i_scale(1) = 48.d0 * dsqrt(10.d0)
         i_scale(2) = 24.d0 * dsqrt(22.d0)
         i_scale(3) = 8.d0 * dsqrt(30.d0)
         i_scale(4) = 16.d0 * dsqrt(165.d0)
         i_scale(5) = 8.d0 * dsqrt(30.d0)
         i_scale(6) = 48.d0 * dsqrt(22.d0)
         i_scale(7) = 48.d0 * dsqrt(10.d0)
         i_scale(8) = 24.d0 * dsqrt(22.d0)
         i_scale(9) = 16.d0 * dsqrt(165.d0)
         i_scale(10) = 48.d0 * dsqrt(1155.d0)
         i_scale(11) = 24.d0 * dsqrt(22.d0)
         i_scale(12) = 48.d0 * dsqrt(55.d0)
         i_scale(13) = 48.d0 * dsqrt(55.d0)
      endif

      istart = 1 
      ierd   = 1

c     write(6,*) ' IFLAGS 71 =', iflags(71) 
c     go to 2000 
      if (iflags(71) .eq. 1) go to 1000 
      if (ispherical .ne. 1) go to 1000

c--------------------------------------------------------------------------
c   Spherical coordinates.
c--------------------------------------------------------------------------

      do ishell = 1, nshells
         iend = istart + nfps(ishell) - 1

c---------------------------------------------------------------------------
c   Angular momentum of the shell determines the number of shell components.
c---------------------------------------------------------------------------

         nsh_coords = 2 * iangular(ishell) + 1
         nctr      = nfps(ishell) / nsh_coords
         if (nsh_coords .eq. 1) then 
            do i = istart, iend
               erd_index(ierd) = ierd
               scalars(ierd)   = 1.d0
               ierd = ierd + 1
            enddo
         else

c----------------------------------------------------------------------------
c   Calculate indices for all elements of the current shell.
c----------------------------------------------------------------------------
         
            if (iangular(ishell) .eq. 1) then    ! p shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (pmap(j)-1)*nctr + i - 1
                     scalars(ierd) = 1.d0
                     ierd = ierd + 1
                  enddo
               enddo
            else if (iangular(ishell) .eq. 2) then    ! d shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (dmap(j)-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do j = 1, nsh_coords
                  do i = 1, nctr
                     scalars(ierd) = d_scale(j)
                     ierd = ierd + 1
                  enddo
               enddo
            else if (iangular(ishell) .eq. 3) then    ! f shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (fmap(j)-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do j = 1, nsh_coords
                  do i = 1, nctr
                     scalars(ierd) = f_scale(j)
                     ierd = ierd + 1
                  enddo
               enddo
            else if (iangular(ishell) .eq. 4) then   ! g shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (gmap(j)-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do j = 1, nsh_coords
                  do i = 1, nctr
                     scalars(ierd) = g_scale(j)
                     ierd = ierd + 1
                  enddo
               enddo
            else if (iangular(ishell) .eq. 5) then   ! h shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (hmap(j)-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do j = 1, nsh_coords
                  do i = 1, nctr
                     scalars(ierd) = h_scale(j)
                     ierd = ierd + 1
                  enddo
               enddo
            else if (iangular(ishell) .eq. 6) then   ! i shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (imap(j)-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do j = 1, nsh_coords
                  do i = 1, nctr
                     scalars(ierd) = i_scale(j)
                     ierd = ierd + 1
                  enddo
               enddo
            else
               print *,'Cannot have shell types > i.'
               call abort_job()
            endif 
         endif

         istart = istart + nfps(ishell)
      enddo
      return

1000  continue

c---------------------------------------------------------------------------
c   Cartesian coordinates.
c---------------------------------------------------------------------------

c --- s-functions --- 
      x_scale(0,1) = 1.0 
c --- p-functions --- 
      do i = 1, 3 
         x_scale(1,i) = 1.0 
      enddo 

c --- d-functions --- 
      x_scale(2,1) = 1.0/dsqrt(3.0d0)  
      x_scale(2,2) = 1.0 
      x_scale(2,3) = 1.0 
      x_scale(2,4) = 1.0/dsqrt(3.0d0)  
      x_scale(2,5) = 1.0 
      x_scale(2,6) = 1.0/dsqrt(3.0d0)  

c --- f-functions --- 
      x_scale(3, 1) = 1.0/dsqrt(15.0d0) 
      x_scale(3, 2) = 1.0/dsqrt(3.0d0) 
      x_scale(3, 3) = 1.0/dsqrt(3.0d0) 
      x_scale(3, 4) = 1.0/dsqrt(3.0d0) 
      x_scale(3, 5) = 1.0/dsqrt(1.0d0) 
      x_scale(3, 6) = 1.0/dsqrt(3.0d0) 
      x_scale(3, 7) = 1.0/dsqrt(15.0d0) 
      x_scale(3, 8) = 1.0/dsqrt(3.0d0) 
      x_scale(3, 9) = 1.0/dsqrt(3.0d0) 
      x_scale(3, 10) = 1.0/dsqrt(15.0d0)

c --- g-functions --- 
      x_scale(4, 1) = 1.0/dsqrt(105.0d0) 
      x_scale(4, 2) = 1.0/dsqrt(15.0d0) 
      x_scale(4, 3) = 1.0/dsqrt(15.0d0) 
      x_scale(4, 4) = 1.0/dsqrt(9.0d0) 
      x_scale(4, 5) = 1.0/dsqrt(3.0d0) 
      x_scale(4, 6) = 1.0/dsqrt(9.0d0) 
      x_scale(4, 7) = 1.0/dsqrt(15.0d0) 
      x_scale(4, 8) = 1.0/dsqrt(8.0d0) 
      x_scale(4, 9) = 1.0/dsqrt(3.0d0) 
      x_scale(4, 10) = 1.0/dsqrt(15.0d0) 
      x_scale(4, 11) = 1.0/dsqrt(105.0d0) 
      x_scale(4, 12) = 1.0/dsqrt(15.0d0) 
      x_scale(4, 13) = 1.0/dsqrt(9.0d0) 
      x_scale(4, 14) = 1.0/dsqrt(15.0d0) 
      x_scale(4, 15) = 1.0/dsqrt(105.0d0) 

c --- h-functions --- 
      x_scale(5, 1) = 1.0/dsqrt(945.0d0) 
      x_scale(5, 2) = 1.0/dsqrt(105.0d0) 
      x_scale(5, 3) = 1.0/dsqrt(105.0d0) 
      x_scale(5, 4) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 5) = 1.0/dsqrt(15.0d0) 
      x_scale(5, 6) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 7) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 8) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 9) = 1.0/dsqrt(9.0d0) 
      x_scale(5, 10) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 11) = 1.0/dsqrt(105.0d0) 
      x_scale(5, 12) = 1.0/dsqrt(15.0d0) 
      x_scale(5, 13) = 1.0/dsqrt(9.0d0) 
      x_scale(5, 14) = 1.0/dsqrt(15.0d0) 
      x_scale(5, 15) = 1.0/dsqrt(105.0d0) 
      x_scale(5, 16) = 1.0/dsqrt(945.0d0) 
      x_scale(5, 17) = 1.0/dsqrt(105.0d0) 
      x_scale(5, 18) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 19) = 1.0/dsqrt(45.0d0) 
      x_scale(5, 20) = 1.0/dsqrt(105.0d0) 
      x_scale(5, 21) = 1.0/dsqrt(945.0d0) 

c --- i-functions --- 
      x_scale(6, 1) = 1.0/dsqrt(10395.0d0) 
      x_scale(6, 2) = 1.0/dsqrt(945.0d0) 
      x_scale(6, 3) = 1.0/dsqrt(945.0d0) 
      x_scale(6, 4) = 1.0/dsqrt(315.0d0) 
      x_scale(6, 5) = 1.0/dsqrt(105.0d0) 
      x_scale(6, 6) = 1.0/dsqrt(315.0d0) 
      x_scale(6, 7) = 1.0/dsqrt(225.0d0) 
      x_scale(6, 8) = 1.0/dsqrt(45.0d0) 
      x_scale(6, 9) = 1.0/dsqrt(45.0d0) 
      x_scale(6, 10) = 1.0/dsqrt(225.0d0) 
      x_scale(6, 11) = 1.0/dsqrt(315.0d0) 
      x_scale(6, 12) = 1.0/dsqrt(45.0d0) 
      x_scale(6, 13) = 1.0/dsqrt(27.0d0) 
      x_scale(6, 14) = 1.0/dsqrt(45.0d0) 
      x_scale(6, 15) = 1.0/dsqrt(315.0d0) 
      x_scale(6, 16) = 1.0/dsqrt(945.0d0) 
      x_scale(6, 17) = 1.0/dsqrt(105.0d0) 
      x_scale(6, 18) = 1.0/dsqrt(45.0d0) 
      x_scale(6, 19) = 1.0/dsqrt(45.0d0) 
      x_scale(6, 20) = 1.0/dsqrt(105.0d0) 
      x_scale(6, 21) = 1.0/dsqrt(945.0d0) 
      x_scale(6, 22) = 1.0/dsqrt(10395.0d0) 
      x_scale(6, 23) = 1.0/dsqrt(945.0d0) 
      x_scale(6, 24) = 1.0/dsqrt(315.0d0) 
      x_scale(6, 25) = 1.0/dsqrt(225.0d0) 
      x_scale(6, 26) = 1.0/dsqrt(315.0d0) 
      x_scale(6, 27) = 1.0/dsqrt(945.0d0) 
      x_scale(6, 28) = 1.0/dsqrt(10395.0d0) 


      do ishell = 1, nshells
         iend = istart + nfps(ishell) - 1

c---------------------------------------------------------------------------
c   Angular momentum of the shell determines the number of shell components.
c---------------------------------------------------------------------------

         nsh_coords = (iangular(ishell)+1)*(iangular(ishell)+2)/2
         nctr      = nfps(ishell) / nsh_coords

         if (iangular(ishell) .lt. 2) then 

         do i = 1, nctr
         do j = 1, nsh_coords
            erd_index(ierd) = istart + (j-1)*nctr + i - 1
            scalars(ierd)   = 1.d0 
            ierd = ierd + 1
         enddo
         enddo

         endif 

         if (iangular(ishell) .ge. 2) then 

         do i = 1, nctr
         do j = 1, nsh_coords
c           xscale = 1.0 
c           if (j .eq. 1 .or. j .eq. 4 .or. j .eq. 6) 
c    #          xscale=1.0/dsqrt(3.0d0) 
            erd_index(ierd) = istart + (j-1)*nctr + i - 1
            scalars(ierd)   = x_scale(iangular(ishell),j)  
            ierd = ierd + 1
         enddo
         enddo

         endif 

c        if (iangular(ishell) .gt. 2) then 

c        do i = 1, nctr
c        do j = 1, nsh_coords
c           xscale = 1.0 
c           erd_index(ierd) = ierd  
c           scalars(ierd)   = 1.d0   
c           ierd = ierd + 1
c        enddo
c        enddo

c        endif 

        istart = istart + nfps(ishell)
      enddo

      return

2000  continue

c---------------------------------------------------------------------------
c   Do not allow any scaling OR reordering  
c   NOT used currently ??? 
c---------------------------------------------------------------------------

      do ishell = 1, nshells
         if (iangular(ishell) .eq. 0) xscale=1.0    ! s shell
         if (iangular(ishell) .eq. 1) xscale=1.0    ! p shell
         if (iangular(ishell) .eq. 2) xscale=1.0/dsqrt(3.0d0) ! d shell
         if (iangular(ishell) .eq. 3) xscale=1.0/dsqrt(4.0d0) ! f shell

c---------------------------------------------------------------------------
c   Angular momentum of the shell determines the number of shell components.
c---------------------------------------------------------------------------

         nsh_coords = (iangular(ishell)+1)*(iangular(ishell)+2)/2
         nctr      = nfps(ishell) / nsh_coords

         do i = 1, nctr
         do j = 1, nsh_coords
            erd_index(ierd) = ierd  
            scalars(ierd)   = xscale  
            write(6,*) ' SCALARS :', ierd, scalars(ierd) 
            ierd = ierd + 1
         enddo
         enddo

      enddo

      return
      end
