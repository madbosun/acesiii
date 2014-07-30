      subroutine rdbasis(
     &     lumol, natm, compnorb, comppopv, compmemb,
     &     nprm, maxshl, nangatm, nshlatm, nshlatm1)
c
c This routine reads VMOL's MOL file to get the following information
c about the atomic basis: NPRM - total number of primitives; MAXSHL -
c the largest number of shells per atom; NSHLATM (NATM)- the number of
c shells for each atom.
c
c These parameters are necessary to allocate memory for the primitive
c exponents and coefficients. After memory is allocated, the loop over
c the MOL file will be repeated to actually read the exponents and
c coefficients.
c
c This routine is a slightly modified copy of LIBRA3/A3BASIS, which in
c turn is a copy of Nevin Oliphant's HFDFT/BASIS routine.
c
c Igor Schweigert, Aug 2003
c $Id: rdbasis.F,v 1.1 2004/12/21 17:41:15 igor Exp $
c
      implicit none
c
c     Arguments.
c     
      integer
     &     lumol, natm, compnorb, compmemb (compnorb), comppopv (natm),
     &  nprm, maxshl, nangatm (natm), nshlatm (natm),nshlatm1 (natm),
     &  ii
c
c     Local variables.
c     
      integer
     &     memb, iorb, ishl, iang, nang, isub, nsub (99),
     &     iprm, iaos, i, nline, ishl1
c     
      character*80 line

c
c     Initialize accumulative quantities to zeros.
c     
      nprm = 0
      maxshl = 0
c
c     Rewind the MOL file and skip the first 5 lines
c     
      rewind(lumol)
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
      read(lumol,'(a)') line
c
c     Loop over symmetry unique atoms. For each unique atom:
c      1/ read number of different angular momenta, and the number of lines
c         describing each angular momentum;
c      2/ loop over angular momentum, and over all lines describing it in order
c         to accumulate the number of shells per atom, and the total number of 
c         primitives (note that the latter is multiplied by the number of
c         symmetry identical atoms corresponding to current symm. unique atom)
c      3/ update the largest number of shells per atom
c      4/ copy the number of shells to to all symmetry identical atoms corresponding
c         to the current symmetry unique atom and shift MEMB by number of symm.
c         identical atoms, COMPPOPV
c
c     Note that in the begging of the loop over symmetry unique atoms
c     (orbits), MEMB points at the first atoms in the block of symmetry
c     identical atoms in COMPMEMB. After the info is copied to all
c     symm. ident. atoms, MEMB is iterated by the current COMPPOPV.
c     
      memb = 1
      do iorb=1, compnorb
c     
         read(lumol, '(a25,12i5)')
     &        line, nang, (nsub(iang), iang=1, nang)
         read(lumol, '(a)') line
c     
         ishl = 0
         ishl1 = 0
         do iang = 1, nang
            ishl = ishl + iang*(iang+1)/2
            ishl1 = ishl1+1
c     
            do isub=1, nsub(iang)
c     
               read(lumol,'(2i5)') iprm, iaos
               nprm = nprm + iprm*iang*(iang+1)/2*comppopv (iorb)
c     
               nline = (iaos-3)/4
               if ((iaos-3).gt.(nline*4)) nline = nline + 1
               nline = (nline+1) * iprm
               do i=1, nline
                  read(lumol,'(a)') line
               end do
c
            end do
         end do
c     
         if (ishl .gt. maxshl) maxshl = ishl
c     
         do i = 1, comppopv (iorb)
            nangatm (compmemb (memb+i-1)) = nang
            nshlatm (compmemb (memb+i-1)) = ishl
            nshlatm1 (compmemb (memb+i-1)) = ishl1
         end do
         memb = memb + comppopv (iorb)
c     
      end do
c 
c      do ii=1,natm
c      write(6,*) nshlatm1 (ii), ii
c      enddo
c      write(6,*) '-----------------------'


    
      return
      end

