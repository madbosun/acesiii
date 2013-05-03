      subroutine compute_return_2SOIX(a1,a2,b1,b2,c1,c2,d1,d2,scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 SOI_JX, nsend)

c---------------------------------------------------------------------------
c   The block of integrals (a1:a2,b1:b2,c1:c2,d1:d2) is computed for the
c   following 'types' of integrals based on atomic labels.
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

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

      INTEGER     DER1X,DER2X,DER3X,DER4X
      INTEGER     DER1Y,DER2Y,DER3Y,DER4Y
      INTEGER     DER1Z,DER2Z,DER3Z,DER4Z

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
      double precision SOI_JX(a1:a2,b1:b2,c1:c2,d1:d2)
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

      common /Imax_com/sz_max(max_shells,max_shells), delta
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

      call mpi_comm_rank(mpi_comm_world, me, ierr)
c      print *,'Task ',me,' computing integrals for ',a1,a2,b1,b2,
c     *     c1,c2,d1,d2
c      call c_flush_stdout()

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1
      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

c Set the integral tolerance

      call set_itol(delta,itol)

c     write(6,*) ' DELTA = ', delta, itol
c     itol = delta ! 1.0d-9

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

CSSS         Write(6,*)
CSSS         Write(6,*) "Loop over shell quads."
CSSS         Write(6,"(4(a,1x,I4,1x,I4,1x))") "m1 and m2: ", m1, m2,
CSSS     &                                    "n1 and N2: ", n1, n2,
CSSS     &                                    "r1 and r2: ", r1, r2,
CSSS     &                                    "s1 and s2: ", s1, s2

C   These are always undifferentiated

CSSS      Write(6,*) "@-entry to OED_batch"
CSSS      Write(6,"((a,2(1x,i4)))") "m, n :", m, n

         DER3X = 0
         DER3Y = 0
         DER3Z = 0
         DER4X = 0
         DER4Y = 0
         DER4Z = 0

c constructing < mu nu | Jz | del gamma >
c more precisely
c Jz = <mu^x nu | del^y gamma> - <mu^y nu | del^x gamma >
c Part A
c <X ..| Y ..>

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
               x4 = coords(1,s)
               y4 = coords(2,s)
               z4 = coords(3,s)
               call pack_coeffs(alphas, ixalpha, coeffs, ixpcoef,
     *                          ncfps, npfps, m, n,
     *                          r, s, alpha_pack, nalpha_pack,
     *                          pcoeff_pack, npcoeff_pack,
     *                          ccbeg, ccend, indx_cc,
     *                          ccbeg_pack, ccend_pack)

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------
       ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)

C
c      Write(6,*) "@-entry to ERD_batch"
c      Write(6,"((a,4(1x,i4)))") "m, n, r, s:", m, n, r, s
c      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
c     &                         nalpha_pack, npcoeff_pack
c      Write(6,"((a,4(1x,i4)))") "Ivanagmom:", Ivangmom(m),
c     &      Ivangmom(n), Ivangmom(r), Ivangmom(s)
c      Write(6,"((a,4(1x,i4)))") "Ncfps:", Ncfps(m),
c     &      Ncfps(n), Ncfps(r), Ncfps(s)
c      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m),
c     &      Npfps(n), Npfps(r), Npfps(s)
c      Write(6,*)
c      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha_pack)
c      Write(6,*)
c      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff_pack)
c      Write(6,*)
c      Write(6,"(4(1x,I4))")(CCbeg_pack(i),i=1,Nalpha_pack)
c      Write(6,*)
c      Write(6,"(4(1x,I4))")(CCend_pack(i),i=1,Nalpha_pack)
c      Write(6,*)
c      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
c     &                       ncsum
c      Write(6,*) "The coordiantes"
c      Write(6,*)
c      Write(6,"(3(1x,F10.5))") x1, y1, z1
c      Write(6,"(3(1x,F10.5))") x2, y2, z2
c      Write(6,"(3(1x,F10.5))") x3, y3, z3
c      Write(6,"(3(1x,F10.5))") x4, y4, z4

        DER1X = 0
        DER1Y = 1
        DER1Z = 0
        DER2X = 0
        DER2Y = 0
        DER2Z = 1

       call ERD__SOI_BATCH_XY(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s),
     *                x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr,
     *                nints,
     *                nfirst, scr)

c---------------------------------------------------------------------------
c   Move the integrals into the output block.
c---------------------------------------------------------------------------

            if (nints .gt. 0) then

c       Write(6,*)"The <munu | lam sig> integral matrix for shell block"
c       Write(6,"(5(1x,F15.10))") (scr(nfirst-1+i),i=1,nints)

               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

               call add_integrals(soi_Jx, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)
            endif

c  <Y ...| X..>

        DER1X = 0
        DER1Y = 0
        DER1Z = 1
        DER2X = 0
        DER2Y = 1
        DER2Z = 0

         call ERD__SOI_BATCH_YX(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr,
     *                nints,
     *                nfirst, scr)

      if (nints .gt. 0) then
               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

      call add_integrals(soi_Jx, a1,a2,b1,b2,c1,c2,d1,d2,
     *                            scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,-1.d0)
      endif

         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m
C
C#ifdef _DEBUG_LVL0
      Write(6,*) "The <munu |SOI_JX| lam sig> integral matrix"
      do d = 1, Nbasis
      Write(6,*)
      do c = 1, Nbasis
      Write(6,*)
      do b = 1, Nbasis
         Write(6, "(5(1x,F15.10))") (SOI_JX(a,b,c,d), a=1, nbasis)
      enddo
      enddo
      enddo
C#endif

      return
      end

