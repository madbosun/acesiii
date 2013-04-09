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
      subroutine contract_ecp_derv_integrals(density,
     *                       da1, da2, db1, db2)
c---------------------------------------------------------------------------
c  Computes the set of OED derivative integrals matching "type" and 
c  and contracts them with density block (da1:da2, db1:db2).
c
c   type = 1 --> NAI integrals
c          2 --> KIN integrals
c          3 --> OVL integrals
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer type
      integer da1, da2, db1, db2
      double precision density(da1:da2, db1:db2)
      double precision z(da1:da2, db1:db2)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif
      integer i, ix, m, n, nfirst, nints, nalpha, ncoeff, iatom, ncsum
      integer component
      integer a, b, a1,a2,b1,b2,lastao, igrad
      integer m1, m2, n1, n2
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer ixderc
      integer ind(2)
      integer flopcount
      logical spherical
      logical dexist 
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2, y
c
c ECP derivative integrals 
c 
      integer natoms 
      integer ndim, j, jj, kk  
      integer Ierror 
      double precision etemp((da2-da1)*(db2-db1))

      common /flags/ iflags
      integer iflags(100)

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      spherical = (ispherical .eq. 1)
      lastao    = end_nfps(nshells)

      ind(1) = 1
      ind(2) = 2

      if (intpkg .ne. flocke_package) then
         print *,'Error: contract_oed_derv_integrals requires ERD ',
     *           'integral package.'
call abort_job()
      endif

      INQUIRE(file='decpint.dat',exist=dExist)

      if (dexist) then

      OPEN (Unit = 88, File = 'decpint.dat', Form = 'Unformatted',
     *      Status = 'Old', Iostat = Ierror)
      rewind 88
#ifdef _DEBUG_LVL0
      write(6,*) ' Reading the ECP derivative integrals '
10    Format(2x,'ATOM NUMBER ',2x,I8)
20    Format(6x,'Coordinate ',2x,I8)
30    Format(12x,'Integrals ',2x,2I8,2x,F16.8)
#endif

      endif 

      call ecp_to_erd(nfps, ivangmom, nshells, ispherical,
     *                         ibuf(ierdind), dbuf(iscale_fac)) 

c--------------------------------------------------------------------------
c   If ECP add ecp derivative integrals  

      read(88) natoms 
      do iatom = 1, natoms 
         do component = 1, 3 
            read(88) ndim 
            read(88) (etemp(i), i=1, ndim*(ndim+1)/2) 

            kk = 0 
            do i = 1, ndim 
            do j = 1, i 
               kk = kk + 1 
               z(i,j) = etemp(kk) ! * scalars(i) * scalars(j)  
               z(j,i) = etemp(kk) ! * scalars(i) * scalars(j)  
            enddo 
            enddo 

            call decp_renorm(ndim, z, ibuf(ierdind), dbuf(iscale_fac))
  
c            do a = a1, a2 
c            do b = b1, b2 
c                z(a,b) = decp_erd(a, b) 
#ifdef _DEBUG_LVL0
c                 Write(6,"(a)") "Adding ECP derivative integrals"
c                 Write(6,"(1x, 2I3, 1x,f12.7)") a, b, decp_erd(iatom,
c    *                                    component, a, b)
#endif
c            enddo 
c            enddo 

c--------------------------------------------------------------------------
c   Contract the density block with the block of integrals.
c--------------------------------------------------------------------------

                 call cwork0222(y,
     *                          density,da1,da2,db1,db2,
     *                          da1,da2,db1,db2, ind,
     *                          z,da1,da2,db1,db2,
     *                          da1,da2,db1,db2,ind(1), ind(2),
     *                          flopcount)

c---------------------------------------------------------------------------
c   Sum into the gradient.
c---------------------------------------------------------------------------

                igrad = (iatom-1)*3 + component 
                gradient_data(igrad) = gradient_data(igrad) + y
         enddo ! icomponent = 1, 3 
      enddo ! iatom = 1, natoms 
      close(88) 

      return
      end

      subroutine decp_renorm(nbasis, z, erd_index, scalars)
c-----------------------------------------------------------------------------
c   Using the shell angular momentum and basis function information, this
c   subroutine calculates an array which maps the ACES integral order to 
c   the ERD integral order.  This array may be used to re-order a block 
c   of integrals calculated by the ERD package into a format
corresponding
c   to the ACES (VMOL-based) integrals.
c
c   Arguments:
c      nfps             Number of basis functions per shell.
c      iangular         Shell type (based on angular momentum) per shell.
c                       (i. e. an S shell = 0, P = 1, etc.)
c      nshells          Number of shells contained in nfps and iangular.
c      ispherical       1 = spherical coordinates, 0 = Cartesian
c      erd_index        An array (output) of nbasis indices.  The ith
c                       value of erd_index is the index in the ERD-based 
c                       system corresponding to index "i" in the ACES-based
c                       system.
c      scalars          An array (output) of scale factors to adjust the
c                       ERD integrals to match VMOL integrals.
c-----------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
c     include 'saved_data.h'
      include 'dropmo.h'
      integer erd_index(*)
      double precision scalars(*)
      double precision z(nbasis,nbasis)

      integer i, j

      integer nbasis, ndim, natoms

      integer map(max_basis)
      logical bexist
      logical dexist
      double precision scale(max_basis)

      do i = 1, nbasis 
      do j = 1, nbasis 
         z(i,j) = z(i,j)*scalars(i)*scalars(j)  
      enddo 
      enddo 

      return 
      end 

