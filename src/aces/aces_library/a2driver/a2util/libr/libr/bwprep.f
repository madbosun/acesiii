
c this procedure checks consistency of bw input and prepares
c renumbering vectors between MO numbers and separate numbering of
c particle and hole orbitals

      subroutine bwprep(nocc,nvrt,iuhf)
      integer iflags,n0,norb,nocc,nvrt,n,iocc1,i,j,k,ir,jr
      dimension nocc(2),nvrt(2),norb(2)
      common /syminf/ nstart,nirrep,irrepy(255,2),dirprd(8,8)
      common /flags/ IFLAGS(100)
      logical bIAmOne





cjp
cjp data for multireference state specific Brillouin-Wigner CC method
cjp coded by Jiri Pittner 1998-2000
cjp
      logical isbwcc,masik,isactive,bwgossip,useeq429,scfrefread
      logical bwwarning
      character*256 bwwarntext
      real*8 ecorrbw0,ecorrbw,epsilon0,fockcontr,denomblow,fockcd
      real*8 heff,heffevalr,heffevali,heffevecl,heffevecr,hdiagcontr
      real*8 fock2elcontr,enerscf,hcore,lambdahomotop,hfakt,diishonset
      real*8 fock2elcontr0,enerscf0,cbwstate,totmaxdenom,heffevecrold
      real*8 intruder,hfaktmax
      integer maxorb,maxref,nref, iref, iocc,iocc0
      integer iphnum,invpnum,invhnum,nbwstates,ibwstate
      integer ibwconvg,internfrom,internto,internnum,internindex
      integer ibwpass
      integer nactive,numactive,ihubaccorr,ihomotop
      integer internfrom1,internto1,internindex1,internnum1
      integer ihefferank, iheffefrom, iheffeto, iheffespin, maxexcit
      integer correctiontype
      integer maxbwwarnings, nproc, myproc
c
      parameter(maxorb=512,maxref=32,maxexcit=9,maxbwwarnings=10)
c      NOTE!!! change of maxorb parameter requires format change
c              and character* change in bwread routine!!!
      parameter(denomblow=1d250)
c

cjp common has been splitted in order to avoid problems
cjp with padding on different 32 and 64 bit architectures

      common/bwccint/isbwcc, masik, nref, iref,iocc(maxorb,maxref,2),
     +     iphnum(maxorb,maxref,2),invpnum(maxorb,maxref,2),
     +     invhnum(maxorb,maxref,2),
     +     isactive(maxorb,2),nbwstates,ibwstate(maxref+1),
     +     internfrom(maxref*(maxref-1)/2,maxref,3),
     +     internto(maxref*(maxref-1)/2,maxref,3),
     +     internindex(maxref*(maxref-1)/2,maxref,3),
     +     internnum(maxref,3),
     +     internfrom1(maxref*(maxref-1)/2,maxref,2),
     +     internto1(maxref*(maxref-1)/2,maxref,2),
     +     internindex1(maxref*(maxref-1)/2,maxref,2),
     +     internnum1(maxref,2),
     +     ibwpass,ibwconvg(maxref),bwgossip,useeq429,
     +     nactive(2),numactive(maxorb,2),ihubaccorr, ihomotop,
     +     iocc0(maxorb,2),scfrefread,
     +     ihefferank(maxref,maxref),iheffefrom(maxexcit,maxref,maxref),
     +     iheffeto(maxexcit,maxref,maxref),
     +     iheffespin(maxexcit,maxref,maxref),
     +     correctiontype,bwwarning(maxbwwarnings),
     +     bwwarntext(maxbwwarnings),nproc,myproc

      common/bwccreal/ecorrbw,epsilon0,cbwstate(maxref+1),
     +     fockcontr(maxorb*(maxorb+1)/2,2),fockcd(maxorb,maxref,2),
     +     heff(maxref,maxref),heffevalr(maxref),heffevali(maxref),
     +     heffevecl(maxref,maxref),heffevecr(maxref,maxref),
     +     hdiagcontr(maxref), fock2elcontr(maxorb,2),
     +     enerscf(maxref), hcore(maxorb,2),
     +     lambdahomotop,hfakt,diishonset,enerscf0,
     +     fock2elcontr0(maxorb,2),ecorrbw0,totmaxdenom,
     +     heffevecrold(maxref,maxref),intruder,hfaktmax

c
c
cjp BRIEF DESCRIPTION OF VARIABLES INTRODUCED FOR THE MR-BWCC ROUTINES
cjp IN FACT, A LOT OF THAT COULD BE USEFUL FOR ANY HILBERT-SPACE MR-CC
c
c
c isbwcc ... flag for doing bwcc calculation
c maxbwwarnings, bwwarning, bwwarntext ... serious warnings will be
c    summarized at the end of xvcc output for the user's' convenience
c ihefferank(jref,iref) ... degree of excitation between jref and iref
c iheffefrom(maxexcit,jref,iref) , iheffeto, iheffespin ... list of
c    indices of that excitation, sorted according to spin and then the
c    indices, numbers stored are defined as effective particle-hole
c    indices of reference iref
c ihomotop ... whether to use homotopic transition to the
c    size-extensivity correction, after which iteration (if .ne.0)
c lambdahomotop scaling factor of the geometrical series of
c    lambda 1->0 transition
c hfakt ... current value of the homotopy parameter
c hfaktmax ... maximal homotopy parameter allowed to consider cc
c    equations converged
c diishonset ... at which value of hfact restart diis convergence acceleration
c masik ... prepare sorted integral file for the program by Masik and stop
c nref ... number of reference configurations
c iref ... current reference configuration and fermi vacuum
c bwgossip ... switch on debugging output
c ibwpass ... routines like newt2 have to be splitted in two passes -
c    construction of Heff and amplitude update after heff is diagonalized
c    for backw. compatibility, instead of introducing a new routine
c    the same routine does different things being called twice with
c    different ibwpass value
c ecorrbw ... correlation energy from BWCC - Heff(iref,iref) ...
c    denominator correction
c ecorrbw0 .... ecorrbw, but not scaled by the homotopic factor hfact
c denomblow ... huge number to cause division underflow - used for
c    zeroing out the internal amplitudes automatically
c nactive(spin): total count of active spinorbitals
c numactive(i=1..nactive,spin): number of i-th active spinorbital
c    in sequential numbering
c isactive(maxorb,spin): belongs given orbital to the active space?
c    for RHF, the beta ones must be initialized to be identical with alpha ones
c iocc(maxorb,1..nref,spin): defines the nref reference configurations
c    for both RHF and UHF:  iocc(i,iref,spin)=0 or 1
c iphnum(orbital no, iref, spin): gives the effective number of orbital
c    (both particle and hole ones are counted starting from 1)
c invpnum(eff.p.orb.,iref,spin): gives true orbital no. from the
c    effective particle one
c invhnum(eff.h.orb.,iref,spin): gives true orbital no. from the
c    effective hole one
c    all these three ones must be in RHF case initialised to be equal in the
c    alpha and beta parts to keep the code unique and simple
c internfrom(sequence counter n,iref,ispin) is for the ab spin case, the
c    other ones have to be iuhf-indexed
c    internfrom, internto: they are first and second index of n-th internal
c    excitation when processing reference given by second index to the array
c
c internindex(sequence counter,iref,ispin) is the position of
c    corresponding denominator in the denominator list
c internnum(iref,ispin) - number of internal excitation in that
c    category = max sequence counter here ispin=1,2,3 for AA,BB,AB
c internfrom1 etc. are analogous quantities for monoexcitations, here
c    ispin =1,2 note for later: all intern.... quantities are irrep-specific!
c fockcontr(findex(i,j),ispin) ... addition to the fock matrix of
c    reference no.1 to obtain the fock matrix of current reference
c    (fermi vacuum)
c fockcd(i,iref,ispin) ... diagonal part of that correction for ref. no. iref.
c hcore(i,ispin) ... one electron diagonal hamiltonian elements
c fock2elcontr(i,ispin) ... 2el contribution to the diagonal fock element
c    used temporarily
c hdiagcontr(iref) ... contribution of differences of HF energies of
c    different Fermi vacua to diagonal Heff elements
c enerscf(iref) ... HF energy of iref-th fermi vacuum
c ihubaccorr ... =1 ... calculate the size extenzivity correction for BWCC
c                =2,3 ... second and third pass of that calculation
c iocc0(maxorb,2) ... like iocc, but for dummy reference configuration
c    corresponding to SCF WF
c fock2elcontr0(maxorb,2) ... like fock2elcontr0 but for dummy reference
c    configuration
c enerscf0 ... like enerscf, but for dummy reference config
c scfrefread ... tells to bwprep that SCF reference has been read from input
c    and should not be generated automatically from nocc(ispin)
c nbwstates ... how many states to average
c ibwstates(1..nbwstates),cbwstates() ... their numbers and coefficients
c correctiontype ... 0=DC,L T2 term is removed/scaled, 1=DC/L term is not
c    removed/scaled
c totmaxdenom ... max 1/denom found for given reference's' fermi vacuum -
c    as indication of possible intruder problem
c intruder ... limit of 1/denom to be considered intruder and its
c    amplitude zeroed
c for parallelization
c nproc ... number of processors (counted from 1)
c myproc ... number of the processor currently executing the code
c    (counted from 1)




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

