
c This routine maps an internal file name (limited to 8 characters)
c like 'JOBARC' to a user-defined (external) filename like '../JOBARC.old'.

cjp - For MR-BWCC, various files are affixed with ".##", in which the hashes
c     represent the 2-digit process id.
cYAU - Our parallel file processing has not been fully designed. Once it is, I
c      imagine this 'feature' will be refined or replaced.

c INPUT
c char*(*) SZINT : the internal file name

c OUTPUT
c char*80 SZEXT   : the external file name for use with OPEN or INQUIRE
c int     ILENGTH : the character length of SZEXT

c#define _DEBUG_GFNAME


      subroutine gfname(szInt,szExt,iLength)
      implicit none

c ARGUMENTS
      character*(*) szInt, szExt
      integer iLength

c EXTERNAL FUNCTIONS
      integer fnblnk
      character*1 achar

c PARAMETERS
      character*1 czSpace
      parameter (czSpace=' ')

c INTERNAL VARIABLES
      character*80 szTmp
      integer iEnd, iEnd1, iTmp, iStat
      logical bExist, bDone

c COMMON BLOCKS


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



c ----------------------------------------------------------------------

      iTmp = 0
c   o assert szInt is not empty
      if (fnblnk(szInt).eq.0) then
         print *, '@GFNAME: Assertion failed.'
         print *, '   szInt = "',szInt,'"'
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)

c ----------------------------------------------------------------------

c   o determine the shortest string length of szInt
      iTmp = min(8,len(szInt))
      iEnd = 0
      do while ((szInt(iEnd+1:iEnd+1).ne.czSpace).and.(iEnd.lt.iTmp))
         iEnd = iEnd + 1
      end do

c   o initialize szExt
      if (iEnd.gt.0) then
         szExt(1:iEnd) = szInt(1:iEnd)
         iLength = iEnd
      else
         iLength = 0
         return
      end if

c   o attempt to map the filename to a user-defined file
      inquire(file='FILES',exist=bExist,err=666,iostat=iStat)
      if (bExist) then
         open(unit=92,file='FILES',
     &        status='OLD',form='FORMATTED',
     &        err=666,iostat=iStat)
         rewind(92,err=666,iostat=iStat)
         bDone = .false.
         do while (.not.bDone)
            read(unit=92,fmt='(a)',
     &           end=200,err=666,iostat=iStat) szTmp
            if ((szTmp(1:iEnd).eq.szInt(1:iEnd)) .and.
     &          (szTmp(iEnd+1:iEnd+1).eq.czSpace)     ) then
               iEnd1 = iEnd+1
               do while ((szTmp(iEnd1+1:iEnd1+1).ne.czSpace).and.
     &                   (iEnd1.lt.80))
                  iEnd1 = iEnd1 + 1
               end do
               iLength = iEnd1-(iEnd+1)
               szExt(1:iLength) = szTmp(iEnd+2:iEnd1)
               bDone = .true.
            end if
         end do
 200     continue
         close(92,status='KEEP',err=666,iostat=iStat)
c     end if (bExist)
      end if

c   o add the process id to the end of various files
      if (.not.masik.and.isbwcc) then
         if (szInt(1:iEnd).eq.'MOINTS'.or.
     &       szInt(1:iEnd).eq.'MOABCD'.or.
     &       szInt(1:iEnd).eq.'FOCKCD'    ) then
            szExt(iLength+1:iLength+1) = '.'
            szExt(iLength+2:iLength+2) = achar(48+iref/10)
            szExt(iLength+3:iLength+3) = achar(48+iref-10*(iref/10))
            iLength = iLength+3
         end if
      end if


      return

c   o FILES I/O error
 666  print *, '@GFNAME: I/O error on FILES'
      print *, '         internal file name = "',szInt,'"'
      print '(/)'
      call aces_io_error('GFNAME',92,iStat)

c     end subroutine gfname
      end

