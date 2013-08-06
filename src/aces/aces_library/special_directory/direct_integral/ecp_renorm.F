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
      subroutine ecp_renorm(nbasis, erd_index, scalars) 
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
c     include 'saved_data.h'
      include 'dropmo.h' 

      integer erd_index(*)
      double precision scalars(*)
      double precision ecp_vmol(max_basis,max_basis)
      double precision etemp(max_basis*max_basis),W,ws

      integer i, j, k, ii, jj, kk  
      integer nctr, nsh_coords

      integer Ierror 

      integer nbasis, ndim, natoms  
      double precision ecp_erd, decp_erd
      common /FECP/ndim, ecp_erd(2500,2500)  
c     common /DECP/natoms, decp_erd(50,3,1000,1000)  

      integer map(max_basis) 
      logical bexist 
      logical dexist 
      double precision scale(max_basis) 

      INQUIRE(file='ecpint.dat',exist=bExist) 

      if (bexist) then 

      OPEN (Unit = 77, File = 'ecpint.dat', Form = 'Unformatted', 
     *      Status = 'Old', Iostat = Ierror) 

      rewind 77 
      read(77) ndim 
      read(77) (etemp(i), i=1, ndim*(ndim+1)/2)  
      close(77) 

      else 
           return 

      endif ! bexist 

#ifdef _DEBUG_LVL0
      write(6,*) ' Renormalizing the ECP integrals '
      write(6,*) ' NDIM = ', ndim 
#endif

      if (ndim .ne. nbasis) then 
         write(6,*) ' The number of functions from ECP read is not 
     *   what is expected:', ndim, nbasis 
      endif 

      do i = 1, ndim 
         ii = erd_index(i) 
         map(ii) = i 
         scale(ii) = scalars(i) 
      enddo 

      do i = 1, ndim  
      do j = 1, i 
           ecp_vmol(i,j) = 0.0  
           ecp_vmol(j,i) = 0.0  
      enddo 
      enddo 

      k = 0 
      do i = 1, ndim  
      do j = 1, i 
           k = k + 1 
           ecp_vmol(i,j) = etemp(k) ! * scalars(ii) * scalars(jj) 
           ecp_vmol(j,i) = etemp(k) ! * scalars(ii) * scalars(jj)   
      enddo 
      enddo 

      do i = 1, ndim  
      do j = 1, ndim  
           ii = map(i) 
           jj = map(j) 
           ecp_erd(i,j) = ecp_vmol(i,j) *scalars(i) * scalars(j)  

#ifdef _DEBUG_LVL0
           write(6,"(a,1x,i2,1x,i2,4(1x,F15.7))")
     &         ' ECP :', i, j, ecp_erd(i,j), ecp_vmol(i,j),
     &                                scalars(i),
     &                                scalars(j)
#endif

      enddo
      enddo

#ifdef _DEBUG_LVL0
      Do i = 1, ndim
         write(6,"(4(2x,F15.8))") (ecp_erd(i,j), j=1,ndim)
         Write(6,*)
         do j = i, ndim
         w  = w +ecp_erd(i, j)
         Ws = ws + ecp_erd(i,j)*ecp_erd(i,j)
         enddo
      enddo
      Write(6,*) "Linear and checksums"
      Write(6, "(2(1x, F20.10))") w, Ws
#endif

      return
      end
