
c this procedure checks consistency of bw input and prepares
c renumbering vectors between MO numbers and separate numbering of
c particle and hole orbitals

      subroutine bwprep(nocc,nvrt,iuhf)
      integer iflags,n0,norb,nocc,nvrt,n,iocc1,i,j,k,ir,jr
      dimension nocc(2),nvrt(2),norb(2)
      common /syminf/ nstart,nirrep,irrepy(255,2),dirprd(8,8)
      common /flags/ IFLAGS(100)
      logical bIAmOne
#include "bwcc.com"

      bIAmOne = .true.

c called after crapsi in order to check bwcc input

      norb(1) = nocc(1)+nvrt(1)
      norb(2) = nocc(2)+nvrt(2)

c check if we are in C1 symmetry
      if (nirrep.gt.1) then
         write(6,*) 'MR-BWCC implemented only for C1 symmetry'
         call aces_exit(1)
      end if

c check array dimension
      if (norb(1).gt.maxorb.or.norb(2).gt.maxorb) then
         write(6,*) 'Too many orbitals, increase MAXORB parameter'
         call aces_exit(1)
      end if

c generate occupation of scf reference if necessary
      if (.not.scfrefread) then
         do iu = 1, 2
            do i = 1, norb(iu)
               iocc0(i,iu) = 0
            end do
         end do
         do iu = 1, 2
            do i = 1, nocc(iu)
               iocc0(i,iu) = 1
            end do
         end do
      end if

c check number of alpha and beta electrons
c here do both spins even for RHF
      do iu = 1, 2
         do j = 1, nref
            n  = 0
            n0 = 0
            do i = 1, norb(iu)
               n  = n  + iocc(i,j,iu)
               n0 = n0 + iocc0(i,iu)
            end do
            if (n.ne.nocc(iu) .or. n0.ne.nocc(iu)) then
               write(6,*)
     &            'inconsistent number of electrons in bwcc input'
               call aces_exit(1)
            end if
         end do
      end do

      do iu = 1, 1+iuhf
c prepare boolean switch isactive
c orbital is active, if there exists at least one pair of reference
c configurations having different occupation of that orbital
         do i = 1, norb(iu)
            isactive(i,iu) = .false.
            iocc1 = iocc(i,1,iu)
            do j = 2, nref
               if (iocc(i,j,iu).ne.iocc1) then
                  isactive(i,iu) = .true.
                  goto 10
               end if
            end do
10          continue
c ISACTIVE MUST be defined for beta orbitals even in the RHF case
c in order to use common rhf/uhf code
            if (iuhf.eq.0) isactive(i,2) = isactive(i,1)
         end do
c prepare nactive, numactive - count and indices of active orbitals
         nactive(iu) = 0
         do i = 1, norb(iu)
            if (isactive(i,iu)) then
               nactive(iu) = nactive(iu)+1
               numactive(nactive(iu),iu) = i
            end if
         end do
      end do

c for RHF copy the beta spin  values
      if (iuhf.eq.0) then
         nactive(2) = nactive(1)
         do i = 1, nactive(1)
            numactive(i,2) = numactive(i,1)
         end do
      end if

c prepare renumbering vectors for later use
      do iu = 1, 1+iuhf
         do j = 1, nref
            ip = 0
            ih = 0
            do i = 1, norb(iu)
               if (iocc(i,j,iu).gt.0) then
                  ih = ih+1
                  iphnum(i,j,iu) = ih
                  invhnum(ih,j,iu) = i
               else
                  ip = ip+1
                  iphnum(i,j,iu) = ip
                  invpnum(ip,j,iu) = i
               end if
               if (iuhf.eq.0) then
                  iphnum(i,j,2)   = iphnum(i,j,1)
                  invhnum(ih,j,2) = invhnum(ih,j,1)
                  invpnum(ip,j,2) = invpnum(ip,j,1)
               end if
            end do
         end do
      end do

c compute mutual excitation level of reference configurations
c and the indices of corresponding internal amplitudes
c zero out effective hamiltonian, so that in the unset elements
c no garbage remains
      do ir = 1, nref
         do jr = 1, nref
            heff(jr,ir) = 0.0
            call calcexcit(ihefferank(jr,ir),iheffefrom(1,jr,ir),
     &                     iheffeto(1,jr,ir),iheffespin(1,jr,ir),
     &                     jr,ir,norb)
c check for identical references!
            if (ir.ne.jr .and. ihefferank(jr,ir).eq.0 ) then
               if (bIAmOne) then
                  write(6,*) '@BWPREP: references no.',ir,jr,
     &                       'are identical!'
               end if
               write(*,*) 'identical references found!'
               call aces_exit(1)
            end if
         end do
      end do

c print out and make a warning if some Heff element will be neglected
      if (bIAmOne) then
         write(6,*)
         write(6,*) 'Excitations between reference configurations'
         write(6,*)
         do ir = 1, nref
            write(6,'(32i2)') (ihefferank(ir,jr),jr=1,nref)
         end do
         write(6,*)
         if (iflags(1).ge.10. or. bwgossip) then
c write out detailed summary
            do ir = 1, nref
            do jr = 1, nref
               write(6,*) 'Heff (',ir,',',jr,'): erank: ',
     &                    ihefferank(ir,jr)
               do i = 1, ihefferank(ir,jr)
                  write(6,*) 'spin ',iheffespin(i,ir,jr),' from ',
     &                       invhnum(iheffefrom(i,ir,jr),jr,
     &                               iheffespin(i,ir,jr)),' to ',
     &                       invpnum(iheffeto(i,ir,jr),jr,
     &                               iheffespin(i,ir,jr))
               end do
               write(6,*)
            end do
            end do
            write(6,*)
            write(6,*)
         end if
      end if

c make check and warning
      do ir = 1, nref
      do jr = 1, nref
         if (ihefferank(jr,ir).gt.2) then
            bwwarning(1) = .true.
cYAU          bwwarntext(1)='@BWPREP-W: There are references mutually more'
cYAU     &       //' than biexcited, some Heff elements will be neglected!'
cYAU          if(bIAmOne)write(6,*) bwwarntext(1)
            if (bIAmOne) then
               write(6,*) '@BWPREP: There are references mutually more',
     &                    ' than biexcited,'
               write(6,*) '           some Heff elements will be',
     &                    ' neglected!'
               write(6,*)
            end if
            goto 123
         end if
      end do
      end do
123   continue

c useeq429=.true. - version to be identically compared to masik's' program
c default should be false ... general implem. of eq.4.28 from Hubac's review
c this is obsolete, this whole switching may be deleted from future version
      useeq429 = .false.

      return
      end

