      subroutine rdprim1(
     &     lumol, natm, nprm, naos, maxshl, compnorb,
     &     comppopv, compmemb, nshlatm, 
     &     nprmshl, iprmshl, naosshl, iaosshl,
     &     alpha, pcoef)
c
c This routine reads the atomic basis primitives' exponents and
c contraction coefficients from VMOL's MOL file. It was written based on
c LIBRA3/A3PRIM which, in turn, is a descendant of Nevin Oliphant's
c HFDFT/PRIM routine.
c
c The routine consist of three block: <> First, we read the MOL file to
c get the numbers of primitive and AOs per each shell on each atom. <>
c Second, we set up the offset arrays that describe where a particular
c shell starts in the big matrices ALPHA and PCOEF. <> Third, we go
c through the MOL file again to read in primitives' exponents and
c contraction coefficients, normalize the orbitals, and store them in
c ALPHA and PCOEF.
c
c The original routine HFDFT/PRIM needed extra scratch arrays to read
c the exponents and coefficients on the first loop over the MOL
c file. The designers of LIBRA3/A3PRIM got rid of these scratch arrays
c by reading the MOL file twice. This also significantly simplifies the
c logic of this routine.
c
c To understand the comment and variables in this routine, it is useful
c to introduce the following notation. A shell is a set of primitives
c (or atomic orbitals) one one atoms that have the same L (angular
c momentum) and Lz (its z-component). All the shells with the same L but
c different Lz's form a 'supershell' (or an angular momentum). Due to
c the spherical symmetry, all the shells in a supershell are identical
c for an atom. Thus the MOL file describes only the supershell. And only
c for symmetry unique atoms. So after we read what we need, we will copy
c the info to all the remaining shells of the supershell and then copy
c it to the all the symmetry identical atoms. That is done inside loop
c over 1..comppopv.
c
c Note also that primitive coefficients will come out with atomic
c orbitals sorted like "s s s s x x x y y y z z z", while all the VMOL
c integrals correspond to AOs sorted like "s s s s x y z x y z x y
c z". So to use the contraction coefficent matrix, one needs to reorder
c it. With some extra effort this can be coded within this procedure (by
c modifying the setup of AOs offset array). But for a time being use the
c PRM2VMOL procedure.
c
c In:  lumol, natm, nprm, naos, maxshl, compnorb, nshlatm, comppopv, compmemb
c Out: nprmshl, iprmshl, naosshl, iaosshl, alpha, pcoef
c
c Igor Schweigert, Sep 2003
c $Id: rdprim.F,v 1.1 2004/12/21 17:41:16 igor Exp $
c
      implicit none
c
c     Input parameters
c     
      integer
     &     lumol, natm, nprm, naos, maxshl, compnorb,
     &     comppopv (natm), compmemb (natm),
     &     nshlatm (natm), 
     &     naosshl (maxshl, natm), nprmshl (maxshl, natm),
     &     iprmshl (maxshl, natm), iaosshl (maxshl, natm)
c     
      double precision
     &     alpha (nprm), pcoef (nprm, naos)
c     
c     Local variables
c     
      integer
     &     memb, iatm, jatm, iorb, iang, nang, isub, nsub (20),
     &     nsubprm, nsubaos, 
     &     ishl, jshl, iprm, iaos, i, j, k
c     
      double precision
     &     x, norm, pi
      parameter (pi = 3.14159265358979d0)
c     
      character*80 line
c     
c     Rewind the MOL file and skip the first 5 lines.
c     
      rewind(lumol)
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
c
c     Read the MOL file to get the number of primitives and AOs per
c     shell. To do this we will loop over the blocks of the symmetry
c     unique atoms (with MEMB pointing at the first atom of the current
c     block). For each block we will: <> Update IATM (points at the
c     first atom in the block) and MEMB (points at the next block).  <>
c     Read in the number of different angular momentum supershells and
c     the number of the MOL file lines describing each supershell. <>
c     Loop over the supershells (with ISHL pointing at the first shell
c     in the current supershell.)  <<>> Loop over the lines describing
c     the current supershell and update the number of primitives and AOs
c     per shell. <> Copy the info to the remaining shells and symmetry
c     identical atoms.
c
      call izero (nprmshl, maxshl*natm)
      call izero (naosshl, maxshl*natm)
c     
      memb = 1
      do iorb=1, compnorb
c
         read(lumol, '(a25,12i5)') line, nang, (nsub(iang), iang=1,nang)
         read(lumol, '(a)') line
c
         iatm = compmemb (memb)
         ishl = 1
         do iang = 1, nang
c
            do isub = 1, nsub(iang)
               read(lumol,'(2i5)') nsubprm, nsubaos
               nprmshl (ishl, iatm) = nprmshl (ishl, iatm) + nsubprm
               naosshl (ishl, iatm) = naosshl (ishl, iatm) + nsubaos
               do i=1, nsubprm
                  read (lumol, '(4f18.10)') x, (x, j=1, nsubaos)
               enddo
            enddo

            do i = 1, comppopv (iorb)
               jatm = compmemb (memb+i-1)
               do jshl = ishl, ishl + iang*(iang+1)/2 - 1
                  nprmshl (jshl, jatm) = nprmshl (ishl, iatm)
                  naosshl (jshl, jatm) = naosshl (ishl, iatm)
               enddo
            enddo
c     
            ishl = ishl + iang*(iang+1)/2
         end do
c     
         memb = memb + comppopv (iorb)
      enddo
c
c     Set up the offset arrays. They will point at where the block of
c     primitives or AOs for a certain shell and atom starts in the big
c     matrix of exponents and primitive coefficients.
c     
      iprm = 1
      iaos = 1
      do iatm = 1, natm
         do ishl = 1, nshlatm (iatm)
            iprmshl (ishl, iatm) = iprm
            iaosshl (ishl, iatm) = iaos
            iprm = iprm + nprmshl (ishl, iatm)
            iaos = iaos + naosshl (ishl, iatm)
         enddo
      enddo
c     
c     Rewind the MOL file again and skip the first 5 lines
c     
      rewind(lumol)
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
c
c     Go through the MOL file again to read in the primitives' exponents
c     and contraction coefficients. To do this we will go over the same
c     loops as in the first time, but now reading in those skipped lines
c     to get exponents and primitive coefficients, and renormalizing the
c     latter ones:
c
      call zero (alpha, nprm)
      call zero (pcoef, nprm*naos)
c     
      memb = 1
      do iorb=1, compnorb
c     
         read(lumol, '(a25,12i5)')
     &        line, nang, (nsub(iang), iang=1, nang)
         read(lumol, '(a)') line
c
         ishl = 1
         do iang = 1, nang
c     
            iprm = iprmshl (ishl, compmemb (memb))
            iaos = iaosshl (ishl, compmemb (memb))
            do isub=1, nsub(iang)
c     
               read(lumol,'(2i5)') nsubprm, nsubaos
c     
               do i=iprm, iprm+nsubprm-1
                  read (lumol, '(4f18.10)')
     &                 alpha (i), (
     &                 pcoef (i,j), j=iaos, iaos+nsubaos-1)
               enddo
c     
               do k=iaos, iaos+nsubaos-1
                  norm = 0.d0
                  do i=iprm, iprm+nsubprm-1
                     do j=1, i
                        x = pcoef (i, k) * pcoef (j, k) *
     &                       (2.d0*dsqrt(alpha (i)*alpha (j)) /
     &                       (alpha (i) + alpha (j)))**
     &                       (dble (iang) + 0.5d0)
                        norm = norm + x
                        if (i .ne. j) norm = norm + x
                     enddo
                  enddo
       
                  do i=iprm, iprm+nsubprm-1
                     pcoef (i, k) = pcoef (i, k) / dsqrt (norm) * 
     &                    (.5d0/pi)**(.75d0) *
     &                    (4.d0*alpha (i))**(.5d0*dble (iang) + .25d0)
                  enddo
               enddo
c     
               iprm = iprm + nsubprm
               iaos = iaos + nsubaos
            enddo
c     
            ishl = ishl + iang*(iang+1)/2
c         write(6,*) iang,"---", dble(iang)
         end do
c     
         do i = 1, comppopv (iorb)
            ishl = 1
            do iang = 1, nang
               do jshl = ishl, ishl+iang*(iang+1)/2-1
c     
                  call dcopy (
     &                 nprmshl (ishl, compmemb (memb)),
     &                 alpha (iprmshl (ishl, compmemb(memb))), 1,
     &                 alpha (iprmshl (jshl, compmemb(memb+i-1))),1)
c     
                  do k=1, naosshl (ishl, compmemb (memb))
                     call dcopy (
     &                    nprmshl (ishl, compmemb (memb)),
     &                    pcoef (
     &                    iprmshl (ishl, compmemb (memb)),
     &                    iaosshl (ishl, compmemb (memb))+k-1), 1,
     &                    pcoef (
     &                    iprmshl (jshl, compmemb (memb+i-1)),
     &                    iaosshl (jshl, compmemb (memb+i-1))+k-1),
     &                    1)
                  enddo
c     
               enddo
               ishl = ishl + iang*(iang+1)/2
            enddo
         enddo
         memb = memb + comppopv (iorb)
c     
      end do
c
      return
      end
