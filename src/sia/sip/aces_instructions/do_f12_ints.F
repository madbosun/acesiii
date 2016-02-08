c---------------------------------------------------------------|
c  This subroutine compute F12-related integrals for given      |
c  block of mu,nu,lambda,sigma superindices                     |     
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine do_f12_ints(a1,a2,b1,b2,c1,c2,d1,d2,scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 out, nsend,int_type)
c---------------------------------------------------------------------------
c   The block of integrals (a1:a2,b1:b2,c1:c2,d1:d2) is computed for the 
c   following 'types' of integrals based on atomic labels.
c---------------------------------------------------------------------------

      implicit none


      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'F12.h'

      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer a,b,c,d

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst
      integer me, ierr
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
      double precision out(a1:a2,b1:b2,c1:c2,d1:d2)
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
      integer int_type

      common /Imax_com/sz_max(max_shells,max_shells), delta
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

#ifdef F12

      call mpi_comm_rank(mpi_comm_world, me, ierr)

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1
      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

      nsend = adim*bdim*cdim*ddim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim,ddim = ',adim,bdim,cdim,ddim
         call mpi_abort(mpi_comm_world, ierr)
      endif

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------


         call lookup_shell(end_nfps, nshells, a1, m1)
         call lookup_shell(end_nfps, nshells, a2, m2)
         call lookup_shell(end_nfps, nshells, b1, n1)
         call lookup_shell(end_nfps, nshells, b2, n2)
         call lookup_shell(end_nfps, nshells, c1, r1)
         call lookup_shell(end_nfps, nshells, c2, r2)
         call lookup_shell(end_nfps, nshells, d1, s1)
         call lookup_shell(end_nfps, nshells, d2, s2)

         do m = m1, m2
            if (m .eq. 1) then
               aa1 = 1
            else
               aa1 = end_nfps(m-1) + 1
            endif
            aa2 = end_nfps(m)

            x1 = coords(1,m)
            y1 = coords(2,m)
            z1 = coords(3,m)
         do n = n1, n2
            if (n .eq. 1) then
               bb1 = 1
            else
               bb1 = end_nfps(n-1) + 1
            endif
            bb2 = end_nfps(n)

            x2 = coords(1,n)
            y2 = coords(2,n)
            z2 = coords(3,n)
         do r = r1, r2
            if (r .eq. 1) then
               cc1 = 1
            else
               cc1 = end_nfps(r-1) + 1
            endif
            cc2 = end_nfps(r)

            x3 = coords(1,r)
            y3 = coords(2,r)
            z3 = coords(3,r)
         do s = s1, s2
            if (s .eq. 1) then
               dd1 = 1
            else
               dd1 = end_nfps(s-1) + 1
            endif
            dd2 = end_nfps(s)

               x4 = coords(1,s)
               y4 = coords(2,s)
               z4 = coords(3,s)
               call pack_coeffs(alphas, ixalpha, coeffs, ixpcoef,
     *                          ncfps, npfps, m, n,
     *                          r, s, alpha_pack, nalpha_pack,
     *                          pcoeff_pack, npcoeff_pack,
     *                          ccbeg, ccend, indx_cc,
     *                          ccbeg_pack, ccend_pack,
     *                          ccbeg_pack64, ccend_pack64)

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------

               ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)



                     IF(int_type.eq.1) then
c--------------------------------------------------------------------------|
c   Here we call cool Norbert's integrator in order to calculate           |
c           Yukawa integrals                                               | 
c--------------------------------------------------------------------------|
c                    write(6,*) 'F12GAMMA=',F12GAMMA
                         

                  call ERD__GENER_YUKAWA_BATCH(intmax, zmax,F12GAMMA,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints,
     *                nfirst, scr)
c---------------------------------------------------------------------------

                      ELSEIF(int_type.eq.2) then
c--------------------------------------------------------------------------|
c   Here we call cool Norbert's integrator in order to calculate           |
c           Slater integrals                                               | 
c--------------------------------------------------------------------------|
c                    write(6,*) 'F12GAMMA=',F12GAMMA
          
                  call ERD__GENER_F12_BATCH(intmax, zmax,F12GAMMA,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints,
     *                nfirst, scr)
c-----------------------------------------------------------------------------

                      ELSEIF(int_type.eq.3) then

c--------------------------------------------------------------------------|
c   Here we call cool Norbert's integrator in order to calculate           |
c           quadratic Slater integrals                                     | 
c--------------------------------------------------------------------------|
c                    write(6,*) 'F12GAMMA=',F12GAMMA

                  call ERD__GENER_F12_BATCH(intmax, zmax,2.d0*F12GAMMA,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints,
     *                nfirst, scr)
c---------------------------------------------------------------------------

                      ELSEIF(int_type.eq.4) then
c--------------------------------------------------------------------------|
c   Here we call cool Norbert's integrator in order to calculate           |
c           [f12,d2/dx2] integrals                                         | 
c--------------------------------------------------------------------------|
                      nints=0
             CALL ERD__GENER_K12_BATCH( intmax,zmax,F12GAMMA,
     &                nalpha_pack, npcoeff_pack, ncsum,
     &                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     &                npfps(m),npfps(n), npfps(r), npfps(s),
     &                ivangmom(m), ivangmom(n),ivangmom(r),ivangmom(s),
     &                x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                2,0,0,
     &                2,0,0,
     &                2,0,0,
     &                2,0,0,
     &                alpha_pack,pcoeff_pack,ccbeg_pack, ccend_pack,
     &                spherical, .true., iscr, nints, nfirst, scr )

c---------------------------------------------------------------------------

                      ELSEIF(int_type.eq.5) then
c--------------------------------------------------------------------------|
c   Here we call cool Norbert's integrator in order to calculate           |
c           [f12,d2/dy2] integrals                                         | 
c--------------------------------------------------------------------------|

                       nints=0
             CALL ERD__GENER_K12_BATCH( intmax,zmax,F12GAMMA,
     &                nalpha_pack, npcoeff_pack, ncsum,
     &                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     &                npfps(m),npfps(n), npfps(r), npfps(s),
     &                ivangmom(m), ivangmom(n),ivangmom(r),ivangmom(s),
     &                x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                0,2,0,
     &                0,2,0,
     &                0,2,0,
     &                0,2,0,
     &                alpha_pack,pcoeff_pack,ccbeg_pack, ccend_pack,
     &                spherical, .true., iscr, nints, nfirst, scr )



c---------------------------------------------------------------------------

                      ELSEIF(int_type.eq.6) then
c--------------------------------------------------------------------------|
c   Here we call cool Norbert's integrator in order to calculate           |
c           [f12,d2/dz2] integrals                                         | 
c--------------------------------------------------------------------------|
                      nints=0
             CALL ERD__GENER_K12_BATCH( intmax,zmax,F12GAMMA,
     &                nalpha_pack, npcoeff_pack, ncsum,
     &                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     &                npfps(m),npfps(n), npfps(r), npfps(s),
     &                ivangmom(m), ivangmom(n),ivangmom(r),ivangmom(s),
     &                x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                0,0,2,
     &                0,0,2,
     &                0,0,2,
     &                0,0,2,
     &                alpha_pack,pcoeff_pack,ccbeg_pack, ccend_pack,
     &                spherical, .true., iscr, nints, nfirst, scr )



c--------------------------------------------------------------------------
                      ENDIF
c                     write(6,*) 'Integral type -',int_type 
c                     write(6,*) 'nints=',nints
c---------------------------------------------------------------------------
c   Move the integrals into the output block.  
c---------------------------------------------------------------------------

            if (nints .gt. 0) then
               

               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

               call move_integrals(out, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst),
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2)
            endif


30       continue
         enddo   ! s
20       continue
         enddo   ! r

10       continue
         enddo   ! n
100      continue
         enddo   ! m

#endif
  
      return
      end



