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
      subroutine compute_atomscf_coeff(watom, scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 nai, kin, ovl,  
     *                 ca,cb) 
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

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
      double precision nai(nc1:nc2,nd1:nd2)
      double precision kin(nc1:nc2,nd1:nd2)
      double precision ovl(nc1:nc2,nd1:nd2)
      double precision H0T(nc1:nc2,nd1:nd2)

      double precision ca(nc1:nc2,nc1:nc2)
      double precision cb(nc1:nc2,nc1:nc2)
      integer imap_p, imap_q, istart, iend
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

      common /Imax_com/sz_max(max_shells,max_shells), delta 
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax    

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter 

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

c      call mpi_comm_rank(mpi_comm_world, me, ierr)

      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

      iatom = watom 
C      write(6,*) ' Performing an SCF calculation on atom:', iatom 

      call comp_return_h0(H0T, iatom, nc1, nc2, nc1, nc2) 

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

         m1 = 1 
         n1 = 1 
         r1 = 1 
         s1 = 1 

         m2 = (nshells)   
         n2 = (nshells)  
         r2 = (nshells) 
         s2 = (nshells)  

c-----------------------------------------------------------------------
c   Find the number of basis functions and shells in the atom.  
c-----------------------------------------------------------------------

         n_basis = 0 
         do m = m1, m2 
            if(atom(m) .eq. iatom) n_basis = n_basis + 
     *                end_nfps(m) - end_nfps(m-1)   
         enddo 
C         write(6,*) ' The number of basis functions on atom', 
C     *                iatom, '=', n_basis  

c-----------------------------------------------------------------------
c   Find the mapping from atom <--> molecule.  
c-----------------------------------------------------------------------

         do n = nc1, nc2 
            map(n) = 0 
            umap(n) = 0 
         enddo 

         n_basis = 0 
         do m = m1, m2 
            beg_anfps(m) = 0   
            end_anfps(m) = 0   

            if(atom(m) .eq. iatom) then 
               beg_anfps(m) = n_basis + 1  

               if (m .eq. 1) then
                  DO n = 1, end_nfps(m) 
                     n_basis = n_basis + 1 
                     map(n_basis) = n 
                     umap(n) = n_basis  
                  enddo 
               else 
                  DO n = end_nfps(m-1) + 1, end_nfps(m) 
                     n_basis = n_basis + 1 
                     map(n_basis) = n 
                     umap(n) = n_basis  
                  enddo 
               endif 

               end_anfps(m) = n_basis 

            endif 

         enddo 

c        do n = nc1, nc2 
c           write(6,*) 'n umap(n)', n, umap(n)  
c        enddo 
c        do m = m1, m2 
c           write(6,*) ' Mth shell:', m, end_anfps(m) 
c        enddo 
C         write(6,*) ' The number of basis functions on atom', 
C     *                iatom, '=', n_basis  

         call do_atomscf_coeff(watom, scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 H0T, nai, kin, ovl,  
     *                 ca, cb,
     *                 n_basis, beg_anfps, end_anfps, map, umap) 

      return 
      end 

      subroutine do_atomscf_coeff(watom, scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 HOT, nai, kin, ovl,  
     *                 cina, cinb,
     *                 n_basis, beg_anfps, end_anfps, map, umap) 
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s, l, mn, rs  
      integer a,b,c,d
      integer iatom, n_basis  
      double precision watom 

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
      double precision nai(nc1:nc2,nd1:nd2)
      double precision kin(nc1:nc2,nd1:nd2)
      double precision ovl(nc1:nc2,nd1:nd2)
      double precision HOT(nc1:nc2,nd1:nd2)
      double precision Cina(nc1:nc2,nd1:nd2)
      double precision Cinb(nc1:nc2,nd1:nd2)

      double precision h0(n_basis,n_basis) 
      double precision aovl(n_basis,n_basis) 
      double precision sos(n_basis,n_basis) 
      double precision Qxx(n_basis,n_basis) 

      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision HFDOLD_A(n_basis,n_basis) 
      double precision HFDOLD_B(n_basis,n_basis) 

      double precision ca(n_basis,n_basis) 
      double precision cb(n_basis,n_basis) 
      double precision cba(n_basis,n_basis) 
      double precision cbb(n_basis,n_basis) 
      double precision FTa(n_basis,n_basis) 
      double precision FTb(n_basis,n_basis) 
      double precision Fa(n_basis,n_basis) 
      double precision Fb(n_basis,n_basis) 
      double precision temp, tempa, tempb   
      integer doit, itemp, p, p1, pp, q, q1, qq

      integer nocc_a, nocc_b, nvirt_a, nvirt_b
      integer iter, max_iter 

      integer imap_p, imap_q, istart, iend
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

      common /Imax_com/sz_max(max_shells,max_shells), delta 
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax    

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter 

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

c     call mpi_comm_rank(mpi_comm_world, me, ierr)

      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

      iatom = watom 

      do m = 1, max_centers 
         if (m .eq. iatom) then 
            nocc_b = charge(m)/2 
            nocc_a = charge(m) - nocc_b 
         endif 
      enddo 

C      write(6,*) ' Performing an SCF calculation on atom:', iatom, 
C     * 'in a basis of', n_basis, 'functions with', nocc_a, nocc_b, 
C     * 'alpha and beta occupied electrons'   

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

      m1 = 1 
      n1 = 1 
      r1 = 1 
      s1 = 1 

      m2 = (nshells)   
      n2 = (nshells)  
      r2 = (nshells) 
      s2 = (nshells)  

c-----------------------------------------------------------------------
c Sum nai and kin into small array and copy ovl there too. 
c --> initial guess   
c-----------------------------------------------------------------------

      itemp = 0 
      do n = nc1, nc2 
      do m = nc1, nc2 
         if (umap(m).ne.0 .and. umap(n).ne.0) then 
            aovl(umap(m),umap(n)) = ovl(m,n) 
c           h0(umap(m),umap(n))   = nai(m,n) + kin(m,n)  
            h0(umap(m),umap(n))   = HOT(m,n)   
            itemp = itemp + 1 
         endif 
      enddo  
      enddo  
      
      if (itemp .ne. n_basis**2) then 
         write(6,*) ' Something wrong with umap ', itemp, n_basis 
         call abort_job()
      endif 

c-----------------------------------------------------------------------
c Construct the hcore initial guess  
c-----------------------------------------------------------------------

      do n = 1, n_basis  
      do m = 1, n_basis  
         FA(m,n) = h0(m,n) 
         FB(m,n) = h0(m,n) 
      enddo  
      enddo  

c-----------------------------------------------------------------------
c Construct U*S**(-1/2)  
c-----------------------------------------------------------------------

      call diag(aovl,sos,m,n_basis,0,1,1) 
      do m = 1, n_basis  
      do n = 1, n_basis  
         temp = 0.0 
         do l = 1, n_basis  
            temp = temp + sos(m,l)*aovl(l,n) 
         enddo 
         Qxx(m,n) = temp 
       enddo  
       enddo  

c-----------------------------------------------------------------------
c Transpose the Fock matrix -> Construct S^(-1/2) F S^(-1/2)  
c-----------------------------------------------------------------------

       call fock_transpose(FA,FB,Qxx,FTa,FTb,n_basis) 

c-----------------------------------------------------------------------
c Diagonalize the transposed Fock matrix  
c-----------------------------------------------------------------------

       call diag(FTa,ca,m,n_basis,0,0,0) 
       call diag(FTb,cb,m,n_basis,0,0,0) 

c-----------------------------------------------------------------------
c Back transform the coefficient array  
c-----------------------------------------------------------------------

       call c_backtran(Qxx,ca,cb,cba,cbb,n_basis) 

c-----------------------------------------------------------------------
c Compute the HF density  
c-----------------------------------------------------------------------

       call hfdensity(ca,cb,HFD_A,HFD_B,n_basis,nocc_a,nocc_b) 

c-----------------------------------------------------------------------
c Compute the HF energy   
c-----------------------------------------------------------------------

       call hfenergy(HFD_A,HFD_B,FA,FB,h0,n_basis) 

c-----------------------------------------------------------------------
c Copy the HF Density into the old Density. 
c-----------------------------------------------------------------------

      call hfdensity_copy(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis)  

c-----------------------------------------------------------------------
c Start the SCF iterations  
c-----------------------------------------------------------------------

      max_iter = 150 !30  
      DO iter = 1, max_iter 

c-----------------------------------------------------------------------
c       Construct the new Fock matrix  
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       One-electron piece  
c-----------------------------------------------------------------------

        do m = 1, n_basis 
        do n = 1, n_basis 
           FA(m,n) = h0(m,n)  
           FB(m,n) = h0(m,n)  
        enddo 
        enddo 

c-----------------------------------------------------------------------
c       Two-electron piece  
c-----------------------------------------------------------------------

         do m = m1, m2
            if(atom(m) .eq. iatom) then  
            aa1 = beg_anfps(m)
            aa2 = end_anfps(m)

            x1 = coords(1,m)
            y1 = coords(2,m)
            z1 = coords(3,m)
         do n = n1, n2
            if (m .le. n) then 
            if(atom(n) .eq. iatom) then  
            bb1 = beg_anfps(n)
            bb2 = end_anfps(n)

            x2 = coords(1,n)
            y2 = coords(2,n)
            z2 = coords(3,n)
         do r = r1, r2
            if(atom(r) .eq. iatom) then  
            cc1 = beg_anfps(r)
            cc2 = end_anfps(r)

            x3 = coords(1,r)
            y3 = coords(2,r)
            z3 = coords(3,r)
         do s = s1, s2
            if (r .le. s) then 
            if(atom(s) .eq. iatom) then  
            dd1 = beg_anfps(s)
            dd2 = end_anfps(s)

            doit = 0 

         if (( r .lt. s) .and.  
     *       ( n .ne. s) .and.  
     *       ( n .ne. r) .and.  
     *       ( m .lt. n) .and.  
     *       ( m .lt. r) .and.  
     *       ( m .ne. s)) doit = 1  

         if ((r .lt. s) .and.  
     *       (m .eq. n)) doit = 1  

         if (( r .lt. s) .and.   
     *    ( n .lt. s) .and.  
     *    ( m .lt. n) .and.  
     *    ( m .eq. r)) doit = 1   

         if (( r .lt. s) .and.  
     *    ( n .eq. r) .and.  
     *    ( m .lt.  n) .and.  
     *    ( m .lt. s)) doit = 1  

         if (( r .lt. s) .and.  
     *    ( n .eq. s) .and.  
     *    ( m .lt. n) .and.  
     *    ( m .lt. r)) doit = 1  
c
         if (( n .eq. m ) .and.  
     *    ( s .eq. m ) .and.  
     *    ( r .eq. m )) doit = 1  
        
         if (( r .eq. s) .and.  
     *    ( m .eq. n) .and.  
     *    ( m .lt. r)) doit = 1  
 
         if (( r .lt. s ) .and.  
     *    ( n .eq. s) .and.  
     *    ( m .lt. n) .and.  
     *    ( m .eq. r)) doit = 1  
 
c
c-----------------------------------------------------------------------
c   Determine the largest density element.
c-----------------------------------------------------------------------

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

c              if (doit .eq. 1) then 

               call ERD__GENER_ERI_BATCH(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum, 
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n), 
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints, 
     *                nfirst, scr)    

c               endif 

c---------------------------------------------------------------------------
c   Move the integrals into the output block.  
c---------------------------------------------------------------------------

           if (nints .gt. 0) then

               call form_ss1fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss2fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss3fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss4fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss5fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss6fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss7fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss8fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

           endif
c 
            endif 
            endif 
         enddo   ! s
            endif 
         enddo   ! r
            endif 
            endif 
         enddo   ! n
            endif 
         enddo   ! m

c-----------------------------------------------------------------------
c       Done computing the Fock matrix   
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       Compute the new HF energy  
c-----------------------------------------------------------------------

        call hfenergy(HFD_A,HFD_B,FA,FB,h0,n_basis) 

c-----------------------------------------------------------------------
c       Transpose the new Fock Matrix   
c-----------------------------------------------------------------------

        call fock_transpose(Fa,FB,Qxx,FTa,FTb,n_basis) 

c-----------------------------------------------------------------------
c       Diagonalize the new Transposed Fock Matrix   
c-----------------------------------------------------------------------

        call diag(FTa,ca,m,n_basis,0,0,0) 
        call diag(FTb,cb,m,n_basis,0,0,0) 

c-----------------------------------------------------------------------
c       Back Transform the coefficient array  
c-----------------------------------------------------------------------

        call c_backtran(Qxx,ca,cb,cba,cbb,n_basis) 

c-----------------------------------------------------------------------
c       Check on convergence and replace the old density with the new   
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       Compute the new HF density  
c-----------------------------------------------------------------------

        call hfdensity(ca,cb,HFD_A,HFD_B,n_basis,nocc_a,nocc_b) 

c-----------------------------------------------------------------------
c       Check for convergence    
c-----------------------------------------------------------------------

        call check_conv(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis,doit)  
        if (doit .eq. 1) go to 100 

c-----------------------------------------------------------------------
c Copy the HF Density into the old Density. 
c-----------------------------------------------------------------------

        call hfdensity_copy(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis)  


      ENDDO ! iter = 1, max_iter 
100   continue 

C      write(6,*) ''
C      write(6,*) ''
C      write(6,*) ' Alpha Orbital energies ' 
C      write(6,*) ' -----------------------------------------------' 
C      do m = 1, nocc_a 
C         write(6,*) '  ', m, FTA(m,m), FTA(m,m)*27.21138386 
C      enddo  
C      write(6,*) ' -----------------------------------------------' 
C      do m = nocc_a + 1, n_basis  
C         write(6,*) '  ', m, FTA(m,m), FTA(m,m)*27.21138386  
C      enddo  
C      write(6,*) ' -----------------------------------------------' 
C      write(6,*) ' ' 
C
C      write(6,*) ' Beta Orbital energies ' 
C      write(6,*) ' -----------------------------------------------' 
C      do m = 1, nocc_b 
C         write(6,*) '  ', m, FTB(m,m), FTB(m,m)*27.21138386 
C      enddo  
C      write(6,*) ' -----------------------------------------------' 
C      do m = nocc_b + 1, n_basis  
C         write(6,*) '  ', m, FTB(m,m), FTB(m,m)*27.21138386  
C      enddo  
C      write(6,*) ' -----------------------------------------------' 
c      write(6,*) ' ' 
c      write(6,*) ' Alpha coefficient matrix '
c      write(6,*) ' -----------------------------------------------' 
c      do m = 1, n_basis
c        write(6,*) '    mu   p   Coefficient '
c        write(6,*) ' -----------------------------------------------' 
c      do n = 1, n_basis
c         if (dabs(ca(n,m)) .gt. 0.0001) then
c         write (6,'(1X,2I5,F12.6)') map(n), map(m), ca (n,m)
c         endif
c      enddo 
c        write(6,*) ' ----------'
c      enddo 
c      write(6,*) ' -----------------------------------------------' 
c      write(6,*) ''
C
C         write(6,*) ' Done SCF calculation of ATOM :', iatom 
C      write(6,*) ''
c      write(6,*) ''

c-----------------------------------------------------------------------
c     Put the coefficient matrix into the full matrix   
c-----------------------------------------------------------------------

   88 FORMAT (A4,I3,A1,I3,A7,I3,A1,I3,A3,F12.6)

      istart = (iatom - 1) * nocc_a + 1
      iend   = istart + nocc_a - 1

      do p = 1, nocc_a
         imap_p = p + istart - 1
         do m = 1, n_basis
            if (map(m) .ne. 0) then
               Cina (map(m),imap_p) = ca (m,p)
C      WRITE (*,88) ' CA(',map(m),',',imap_p,') = ca(',m,',',p,') =',
C     +              ca (m,p)
            endif
         enddo
C         write (*,*) ''
      enddo

      nvirt_a = n_basis - nocc_a
      istart = nalpha_occupied + ((iatom - 1) * nvirt_a + 1)

      do p = 1, nvirt_a
         imap_p = p + istart - 1
         pp = p + nocc_a
         do m = 1, n_basis
            if (map(m) .ne. 0) then
               Cina (map(m),imap_p) = ca (m,pp)
C      WRITE (*,88) ' CA(',map(m),',',imap_p,') = ca(',m,',',pp,') =',
C     +              ca (m,pp)
            endif
         enddo
C         write (*,*) ''
      enddo

      istart = (iatom - 1) * nocc_b + 1
      iend   = istart + nocc_b - 1

      do q = 1, nocc_b
         imap_q = q + istart - 1
         do m = 1, n_basis
            if (map(m) .ne. 0) then
               Cinb (map(m),imap_q) = cb (m,q)
C      WRITE (*,88) ' CB(',map(m),',',imap_q,') = cb(',m,',',q,') =',
C     +              cb (m,q)
            endif
         enddo
C         write (*,*) ''
      enddo

      nvirt_b = n_basis - nocc_b
      istart = nbeta_occupied + ((iatom - 1) * nvirt_b + 1)

      do q = 1, nvirt_b
         imap_q = q + istart - 1
         qq = q + nocc_b
         do m = 1, n_basis
            if (map(m) .ne. 0) then
               Cinb (map(m),imap_q) = cb (m,qq)
C      WRITE (*,88) ' CB(',map(m),',',imap_q,') = cb(',m,',',qq,') =',
C     +              cb (m,qq)
            endif
         enddo
C         write (*,*) ''
      enddo

      return
      end
 
