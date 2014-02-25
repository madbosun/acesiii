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
      subroutine f12_gen_1e(z, da1, da2, db1, db2)
c-----------------------------------------------------------------------------|
c  Computes the set of OED integrals for calculation of 3-electron integrals  |
c-----------------------------------------------------------------------------|
      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'F12.h'
      include 'F12_1el_type.h'
      include 'F12_gr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nsend 
      integer da1, da2, db1, db2
      integer adim, bdim 
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
      integer ierr 

      integer flopcount,iqq
      logical spherical
      logical*8 l8true, l8spherical
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers),CHARGE1(1)
      double precision x1,y1,z1,x2,y2,z2, y,GridX,GridY,GridZ

      double precision RYSDATA(1000000),RYSDATA1(1000000)

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter, dcoord

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif
      spherical = (ispherical .eq. 1)
      l8spherical = spherical
      l8true = .false. ! .true.
      lastao = end_nfps(nshells)

c --------------------------------------------------------------------------- 
c  Set up send data. 
c --------------------------------------------------------------------------- 

      adim = da2-da1+1
      bdim = db2-db1+1 

      nsend = adim*bdim
      if (nsend .lt. 0) then
         print *,'ERROR IN COMP_FOCK_DER ',me,' nsend = ',nsend
         print *,'adim,bdim = ',adim,bdim
         call abort_job()
      endif 

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do b = db1,db2
      do a = da1,da2
         Z(a,b) = 0.d0
      enddo
      enddo

      ind(1) = 1
      ind(2) = 2

      if (intpkg .ne. flocke_package) then
         print *,'Error: comp_fock_der requires ERD ',
     *           'integral package.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Find the shells we need to loop over for this block.
c--------------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, da1, m1)
      call lookup_shell(end_nfps, nshells, da2, m2)
      call lookup_shell(end_nfps, nshells, db1, n1)
      call lookup_shell(end_nfps, nshells, db2, n2)
      
c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = 1, nshells
         iatom = atom(m)
         xn(iatom) = dbuf(icoord+3*(m-1))
         yn(iatom) = dbuf(icoord+3*(m-1)+1)
         zn(iatom) = dbuf(icoord+3*(m-1)+2)
      enddo

c-------------------------------------------------------------------------
c   Calculate the integrals over the necessary shell blocks.
c-------------------------------------------------------------------------

       do m = m1, m2
       do n = n1, n2

            ncsum = ncfps(m) + ncfps(n) 
            x1 = dbuf(icoord+3*(m-1))
            y1 = dbuf(icoord+3*(m-1)+1)
            z1 = dbuf(icoord+3*(m-1)+2)
 
            x2 = dbuf(icoord+3*(n-1))
            y2 = dbuf(icoord+3*(n-1)+1)
            z2 = dbuf(icoord+3*(n-1)+2)

            call pack_coeffs_oed(dbuf(ialpha), ixalpha, dbuf(ipcoeff), 
     *                       ixpcoef, ncfps, npfps, m, n,
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ibuf(iccbeg), ibuf(iccend), 
     *                       indx_cc,
     *                       ccbeg_pack, ccend_pack,
     *                       max_dim_coeff)

c ------------------------------------------------------------------------ 
                            CHARGE1(1) = -1.d0

c---------------Definition of the coordinates of current grid point-------
                   GridX = PTX(cnt)
                   GridY = PTY(cnt)
                   GridZ = PTZ(cnt)

c          write(6,*) '==>',cnt, type_1el,'<==='
c ------------------------------------------------------------------------ 
c        Calculation of <chi_{i}| 1/|r1-rg| |chi_{j}> integrals          |
c ------------------------------------------------------------------------ 
                 IF(type_1el.eq.1) then
                      
             call oed__gener_nai_batch_new(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 1, GridX,GridY,GridZ, CHARGE1, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))
                 ENDIF

c------------------------------------------------------------------------ 
c        Calculation of <chi_{i}| f_{1g} |chi_{j}> integrals            |
c------------------------------------------------------------------------ 


                 IF(type_1el.eq.2) then

             call oed__gener_nai_batch_f12(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 1, GridX,GridY,GridZ, CHARGE1, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), 
     *                 F12GAMMA,1,1,RYSDATA,RYSDATA1,
     *                 nints, nfirst,
     *                 dbuf(dpkgscr))
                 ENDIF

c------------------------------------------------------------------------ 
c        Calculation of <chi_{i}| f_{1g}^2 |chi_{j}> integrals          |
c------------------------------------------------------------------------ 

                 IF(type_1el.eq.3) then

             call oed__gener_nai_batch_f12(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 1, GridX,GridY,GridZ, CHARGE1, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr),
     *                 2.d0*F12GAMMA,1,1,RYSDATA,RYSDATA1,
     *                 nints, nfirst,
     *                 dbuf(dpkgscr))
                 ENDIF

c--------------------------------------------------------------------------|
c      Calculation of <chi_{i}| d(f_{1g})/dx |chi_{j}> integrals           |
c--------------------------------------------------------------------------|

                 IF(type_1el.eq.4) then
             call oed__gener_nai_derv_f12_batch(intmax, zmax, nalpha,
     *            ncoeff,ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *            ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *            1, GridX,GridY,GridZ, CHARGE1,1,
     *            0,0,0,
     *            0,0,0,
     *            1,0,0,
     *            alpha_pack,pcoeff_pack, ccbeg_pack, ccend_pack,
     *            spherical, .true., ibuf(ipkgscr),F12GAMMA,1,1,
     *            RYSDATA,RYSDATA1,
     *            nints, nfirst,
     *            dbuf(dpkgscr))

                 ENDIF

c--------------------------------------------------------------------------|
c      Calculation of <chi_{i}| d(f_{1g})/dx |chi_{j}> integrals           |
c--------------------------------------------------------------------------|

                 IF(type_1el.eq.5) then
             call oed__gener_nai_derv_f12_batch(intmax, zmax, nalpha,
     *            ncoeff,ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *            ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *            1, GridX,GridY,GridZ, CHARGE1,1,
     *            0,0,0,
     *            0,0,0,
     *            0,1,0,
     *            alpha_pack,pcoeff_pack, ccbeg_pack, ccend_pack,
     *            spherical, .true., ibuf(ipkgscr),F12GAMMA,1,1,
     *            RYSDATA,RYSDATA1,
     *            nints, nfirst,
     *            dbuf(dpkgscr))

                 ENDIF


c--------------------------------------------------------------------------|
c      Calculation of <chi_{i}| d(f_{1g})/dx |chi_{j}> integrals           |
c--------------------------------------------------------------------------|

                 IF(type_1el.eq.6) then
             call oed__gener_nai_derv_f12_batch(intmax, zmax, nalpha,
     *            ncoeff,ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *            ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *            1, GridX,GridY,GridZ, CHARGE1,1,
     *            0,0,0,
     *            0,0,0,
     *            0,0,1,
     *            alpha_pack,pcoeff_pack, ccbeg_pack, ccend_pack,
     *            spherical, .true., ibuf(ipkgscr),F12GAMMA,1,1,
     *            RYSDATA,RYSDATA1,
     *            nints, nfirst,
     *            dbuf(dpkgscr))

                 ENDIF

c--------------------------------------------------------------------------| 
c     Calculation of <[d/dx chi_{i}] | d(f_{1g})/dx |chi_{j}> integrals    |          
c--------------------X-derivative------------------------------------------|

                 IF(type_1el.eq.7) then
             call oed__gener_nai_derv_f12g_batch(intmax, zmax, nalpha, 
     *            ncoeff,ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *            ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *            1, GridX,GridY,GridZ, CHARGE1,1, 
     *            1,0,0,
     *            0,0,0,
     *            1,0,0,
     *            alpha_pack,pcoeff_pack, ccbeg_pack, ccend_pack,
     *            spherical, .true., ibuf(ipkgscr),F12GAMMA,1,1,
     *            RYSDATA,RYSDATA1,
     *            nints, nfirst,
     *            dbuf(dpkgscr))
                 ENDIF

c--------------------------------------------------------------------------| 
c     Calculation of <[d/dy chi_{i}] | d(f_{1g})/dy |chi_{j}> integrals    |          
c--------------------Y-derivative------------------------------------------|

                 IF(type_1el.eq.8) then
             call oed__gener_nai_derv_f12g_batch(intmax, zmax, nalpha, 
     *            ncoeff,ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *            ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *            1, GridX,GridY,GridZ, CHARGE1,1, 
     *            0,1,0,
     *            0,0,0,
     *            0,1,0,
     *            alpha_pack,pcoeff_pack, ccbeg_pack, ccend_pack,
     *            spherical, .true., ibuf(ipkgscr),F12GAMMA,1,1,
     *            RYSDATA,RYSDATA1,
     *            nints, nfirst,
     *            dbuf(dpkgscr))
                 ENDIF


c--------------------------------------------------------------------------| 
c     Calculation of <[d/dz chi_{i}] | d(f_{1g})/dz |chi_{j}> integrals    |          
c--------------------Z-derivative------------------------------------------|

                 IF(type_1el.eq.9) then
             call oed__gener_nai_derv_f12g_batch(intmax, zmax, nalpha, 
     *            ncoeff,ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *            ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *            1, GridX,GridY,GridZ, CHARGE1,1, 
     *            0,0,1,
     *            0,0,0,
     *            0,0,1,
     *            alpha_pack,pcoeff_pack, ccbeg_pack, ccend_pack,
     *            spherical, .true., ibuf(ipkgscr),F12GAMMA,1,1,
     *            RYSDATA,RYSDATA1,
     *            nints, nfirst,
     *            dbuf(dpkgscr))
                 ENDIF

c---------------------------------------------------------------------------
               if (nints .gt. 0) then
                  a2 = end_nfps(m)
                  if (m .eq. 1) then
                     a1 = 1
                  else
                     a1 = end_nfps(m-1)+1
                  endif

                  b2 = end_nfps(n)
                  if (n .eq. 1) then
                     b1 = 1
                  else
                     b1 = end_nfps(n-1)+1
                  endif

                 call add_integrals2(z, da1,da2,db1,db2,
     *                 dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2,1.0d0)

               endif ! nints 

         enddo   ! n shells
         enddo   ! m shells

      return
      end
