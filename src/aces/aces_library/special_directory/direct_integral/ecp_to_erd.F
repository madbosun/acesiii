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
      subroutine ecp_to_erd(nfps, iangular, nshells, ispherical,
     *                      erd_index, scalars) 
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
      double precision d_scale(5)
      double precision f_scale(7)
      double precision g_scale(9)
      double precision h_scale(11)
      double precision i_scale(13)

      integer istart, iend, ishell
      integer i, j, ierd
      integer nctr, nsh_coords
  
      double precision x_scale(0:6,1:28) 

      common /flags/ iflags
      integer iflags(100)

      Istart = 1
      Ierd   = 1
     
      If (iSpherical .Eq. 1) Then

         d_scale(1) = 1.0d0 / 2.0D0
         d_scale(2) = 1.d0
         d_scale(3) = 1.0d0 /Dsqrt(3.0d0)
         d_scale(4) = 1.0D0
         d_scale(5) = 1.0D0 / 2.0D0

         f_scale(1) = 1.0d0 / Dsqrt(24.0d0)
         f_scale(2) = 1.0d0 / dsqrt(4.0d0)
         f_scale(3) = 1.0d0 / dsqrt(10.0d0)
         f_scale(4) = dsqrt(6.0d0) / dsqrt(10.0d0)
         f_scale(5) = 1.0d0 / dsqrt(10.0d0)
         f_scale(6) = 1.0d0 / dsqrt(4.0d0)
         f_scale(7) = 1.0d0 / dsqrt(24.0d0)

         g_scale(1) = 1.0d0 / dsqrt(192.0d0)
         g_scale(2) = 1.0d0 / dsqrt(24.0d0)
         g_scale(3) = 1.0d0 / dsqrt(84.0d0)
         g_scale(4) = Dsqrt(18.0d0) / dsqrt(84.0d0)
         g_scale(5) = Dsqrt(3.0d0) / dsqrt(35.0D0)
         g_scale(6) = Dsqrt(18.0d0) / dsqrt(84.0d0)
         g_scale(7) = 1.0d0 / dsqrt(84.0d0)
         g_scale(8) = 1.0d0 / dsqrt(24.0d0)
         g_scale(9) = 1.0d0 / dsqrt(192.0d0)

         h_scale(1) = 1.0D0 / dsqrt(1920.0d0)
         h_scale(2) = 1.0d0 / dsqrt(192.0d0)
         h_scale(3) = 1.0d0 / dsqrt(864.0d0)
         h_scale(4) = 1.0d0 / dsqrt(36.0d0)
         h_scale(5) = 1.d0 / dsqrt(63.0d0)
         h_scale(6) = 1.d0 / dsqrt(4.20d0)
         h_scale(7) = 1.0d0 / dsqrt(63.0d0)
         h_scale(8) = 1.0d0 / dsqrt(36.d0)
         h_scale(9) = 1.0d0 / dsqrt(864.0d0)
         h_scale(10) = 1.0d0 / dsqrt(192.0d0)
         h_scale(11) = 1.0D0 / dsqrt(1920.0d0)

         i_scale(1) = 1.0d0 / dsqrt(23040.0d0)
         i_scale(2) = 1.0d0 / dsqrt(1920.0d0)
         i_scale(3) = 1.0d0 / dsqrt(10560.0d0)
         i_scale(4) = 1.0d0 / dsqrt(352.0d0)
         i_scale(5) = 1.0d0 / dsqrt(792.0d0)
         i_scale(6) = 1.0d0 / dsqrt(19.80d0)
         i_scale(7) = 1.0d0 / dsqrt(46.20d0)
         i_scale(8) =  1.0d0 / dsqrt(19.80d0)
         i_scale(9) =  1.0d0 / dsqrt(792.0d0)
         i_scale(10) = 1.0d0 / dsqrt(352.0d0)
         i_scale(11) = 1.0d0 / dsqrt(10560.0d0)
         i_scale(12) = 1.0d0 / dsqrt(1920.0d0)
         i_scale(13) = 1.0d0 / dsqrt(23040.0d0)
C
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
     *                                 (j-1)*nctr + i - 1
                     scalars(ierd) = 1.d0
                     ierd = ierd + 1
                  enddo
               enddo

            else if (iangular(ishell) .eq. 2) then    ! d shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (j-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do i = 1, nctr
                  do j = 1, nsh_coords
                     scalars(ierd) =  d_scale(j)
                     ierd = ierd + 1
                  enddo
               enddo

            else if (iangular(ishell) .eq. 3) then    ! f shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (j-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do i = 1, nctr
                  do j = 1, nsh_coords
                     scalars(ierd) = f_scale(j)

                     ierd = ierd + 1
                  enddo
               enddo

            else if (iangular(ishell) .eq. 4) then   ! g shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (j-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do i = 1, nctr
                  do j = 1, nsh_coords
                     scalars(ierd) = g_scale(j)

                     ierd = ierd + 1
                  enddo
               enddo

            else if (iangular(ishell) .eq. 5) then   ! h shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (j-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart
               do i = 1, nctr
                  do j = 1, nsh_coords
                     scalars(ierd) = h_scale(j)

                     ierd = ierd + 1
                  enddo
               enddo

            else if (iangular(ishell) .eq. 6) then   ! i shell
               do i = 1, nctr
                  do j = 1, nsh_coords
                     erd_index(ierd) = istart + 
     *                                 (j-1)*nctr + i - 1
                     ierd = ierd + 1
                  enddo
               enddo

               ierd = istart

               do i = 1, nctr
                  do j = 1, nsh_coords
                     scalars(ierd) = i_scale(j)

                     ierd = ierd + 1
                  enddo
               enddo
            else
               print *,'Cannot have shell types > i.'
               call abort_job()
            endif 
C
            Endif

            istart = istart + nfps(ishell)
         enddo

         return

      Else
C
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
                  erd_index(ierd) = istart + (j-1)*nctr + i - 1
                  scalars(ierd)   = x_scale(iangular(ishell),j)  
                  ierd = ierd + 1
               enddo
               enddo

            endif 

            istart = istart + nfps(ishell)
 
         enddo

         return

      Endif

      return
      end
