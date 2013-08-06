
c This routine finalizes the I/O subsystem.

      subroutine aces_io_fin
      implicit none

c EXTERNAL FUNCTIONS
      integer iiamax

c INTERNAL VARIABLES
      integer iUnit, iLenMOIO, iStat
      logical bOpened

c COMMON BLOCKS
c lists.com : begin

c These common blocks contain global information about the arrays in storage.
c Elements prepended with "bw" are for storing the file metadata while working
c on multiple references.






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



      external aces_bd_lists

c moio  (iGrp,iFam) : the physical record that contains the first element
c                     of the array (iGrp,iFam)
c moiowd(iGrp,iFam) : the integer-word index of the first element
c moiods(iGrp,iFam) : the number of columns in the array
c moiosz(iGrp,iFam) : the number of rows    in the array
c moiofl(iGrp,iFam) : the external file unit that contains the array

      integer          moio  (10,500),
     &                 moiowd(10,500),
     &                 moiosz(10,500),
     &                 moiods(10,500),
     &                 moiofl(10,500),
     &               bwmoio  (10,500,maxref),
     &               bwmoiowd(10,500,maxref),
     &               bwmoiosz(10,500,maxref),
     &               bwmoiods(10,500,maxref),
     &               bwmoiofl(10,500,maxref)
      common /lists/   moio,   moiowd,   moiosz,   moiods,   moiofl,
     &               bwmoio, bwmoiowd, bwmoiosz, bwmoiods, bwmoiofl
      save   /lists/

c moiomxsz(iGrp,iFam) : the original length of a one-dimensional array
c                       (This is shameful. Arrays should not be re-dimensioned
c                        at will during a job.)

      integer               moiomxsz(10,500),
     &                    bwmoiomxsz(10,500,maxref)
      common /lists_mxsz/   moiomxsz,
     &                    bwmoiomxsz
      save   /lists_mxsz/

c pRec(i)    : the index of the physical record in file i containing free space
c              (i is the internal unit number of the storage file.)
c iIntOff(i) : the integer offset from the beginning of the physical record
c              needed to address the free space

      integer            pRec   (5),
     &                   iIntOff(5),
     &                 bwpRec   (5,maxref),
     &                 bwiIntOff(5,maxref)
      common /io_ptrs/   pRec,   iIntOff,
     &                 bwpRec, bwiIntOff
      save   /io_ptrs/

c bIOUp  : a flag for bombing in get/putlst if aces_io_init has not been called
c bIOMod : a flag for updating the records in aces_io_fin

      logical           bIOUp, bIOMod
      common /io_flags/ bIOUp, bIOMod
      save   /io_flags/

c lists.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c cache.com : begin

c These common blocks contain global information about the automatic file cache.
c getlst and putlst REQUIRE a cache, hence the term 'automatic' (compared to the
c auxiliary cache controlled by /auxcache/quikget).

c#define _CACHE_BYPASS /* bypasses cache on reading/writing of full records */
c#define _CACHE_HIST
c#define _CACHE_HIST_VERBOSE

      external aces_bd_cache

c icache     : the anchor used to address each cache slot
c cachnum    : the number of usable cache slots
c cachrec(i) : the index of the physical record cached by the data in slot i
c cachfil(i) : the external file unit number that stores the data in slot i
c cachndx(i) : the icache index of slot i
c cachmod(i) : a modification flag used to trigger a writeback

c cachetime   : a cache-event counter
c lrustats(i) : the last 'time' slot i was accessed

      integer        icache(1), cachnum,
     &               cachrec(128),
     &               cachfil(128),
     &               cachndx(128),
     &               cachmod(128)
      common /cache/ icache, cachnum,
     &               cachrec,
     &               cachfil,
     &               cachndx,
     &               cachmod
      save   /cache/

      integer           cachetime, lrustats(128)
      common /cachelru/ cachetime, lrustats
      save   /cachelru/

c cachemiss      : measures cache misses
c cacheskip      : measures read and write bypasses
c cacheread      : measures read hits
c cachewrite     : measures write hits
c cachewriteback : measures writes-back of dirty slots

      integer             cachemiss, cacheskip,
     &                    cacheread, cachewrite, cachewriteback
      common /cache_hist/ cachemiss, cacheskip,
     &                    cacheread, cachewrite, cachewriteback
      save   /cache_hist/

c bCacheUp : a flag for bombing in get/putlst if there is no I/O cache

      logical              bCacheUp
      common /cache_flags/ bCacheUp
      save   /cache_flags/

c cache.com : end
c auxcache.com : begin

c The auxiliary cache is a programmer-controlled list cache. The dimensions are
c the same as those of the MOIO arrays. The programmer may load lists into icore
c memory and set quikget(?,?) to their icore addresses. When getlst and putlst
c operate on ANY list, quikget is checked to see if the list lives in icore.
c If so, the operation is performed on the in-core data instead of hitting the
c storage file(s).

c WARNING
c    There is no automatic updating of storage files with data in the auxiliary
c cache. If the memory-resident data is altered and must be stored to disk, then
c the quikget value must be destroyed (zeroed) and the actual location of the
c data must then be passed to putlst. The aces_auxcache_flush routine
c systematically stores the quikget values, calls putlst, and restores the
c quikget values.


      integer           quikget(10,500)
      common /auxcache/ quikget
      save   /auxcache/

c auxcache.com : end

c ----------------------------------------------------------------------

      iUnit = 0
c   o assert quikget is empty
      iLenMOIO = iiamax(10*500,quikget,1)
      if ((iLenMOIO.ne.1).or.(quikget(1,1).ne.0)) then
         print *, '@ACES_IO_FIN: Assertion failed.'
         print *, '   quikget col = ',1+((iLenMOIO-1)/10)
         print *, '   quikget row = ',1+mod((iLenMOIO-1),10)
         iUnit = 1
      end if
      if (iUnit.ne.0) call aces_exit(iUnit)

c ----------------------------------------------------------------------

c   o make sure the I/O subsystem is up
      if (.not.bIOUp) return

c   o flush the automatic cache
      call aces_cache_flush


c   o update storage records only if a list has been modified
      if (bIOMod) then

c      o get the total number of arrays
         iLenMOIO = 10 * 500

c      o finalize the lists common block
         call putrec(1,'JOBARC','MOIOVEC',iLenMOIO,moio)
         call putrec(1,'JOBARC','MOIOWRD',iLenMOIO,moiowd)
         call putrec(1,'JOBARC','MOIOSIZ',iLenMOIO,moiosz)
         call putrec(1,'JOBARC','MOIODIS',iLenMOIO,moiods)
         call putrec(1,'JOBARC','MOIOFIL',iLenMOIO,moiofl)
         call putrec(1,'JOBARC','MOIOMXSZ',iLenMOIO,moiomxsz)

c      o write out the distribution types
         call putrec(1,'JOBARC','ISYMTYP',2*500,isytyp)

c      o finalize the io_ptrs common block
         call putrec(1,'JOBARC','TOTRECMO',5,pRec)
         call putrec(1,'JOBARC','TOTWRDMO',5,iIntOff)

c      o reset the I/O modification flag
         bIOMod = .false.

c     end if (bIOMod)
      end if


c   o close the external file units
      do iUnit = 50, 50-1+5
         inquire(unit=iUnit,opened=bOpened,err=666,iostat=iStat)
         if (bOpened) then
            close(unit=iUnit,status='KEEP',err=666,iostat=iStat)
         end if
      end do

c   o turn off the I/O subsystem flag
      bIOUp = .false.

      return

c   o I/O error
 666  print *, '@ACES_IO_FIN: I/O error'
      print '(/)'
      call aces_io_error('ACES_IO_FIN',iUnit,iStat)

c     end subroutine aces_io_fin
      end

