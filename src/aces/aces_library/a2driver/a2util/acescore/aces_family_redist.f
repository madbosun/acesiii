
c This routine redistributes the array sizes for an MOIO family of arrays
c based on two-particle symmetry processing. It DOES NOT perform the resort,
c and it requires the groups to be 'packed' in the sense that aces_list_touch
c was not told to start each group on a physical record boundary (besides the
c first one). Although it updates the auxiliary cache addresses, every
c array in the family MUST be in memory with no spaces between them (just
c like the disk storage requirement).
c    This is really a shameful routine. The main idea is that a diagram
c like <Ab|Ij> is to be resorted and stored as <AI|bj>. The programmer
c does the resort with sstgen, calls this routine, and then calls putlst
c to update the disk records. People should consider better algorithms...

c INPUT
c int IRP_DIAG : the total symmetry of the diagram (normally 1)
c int IFAMNDX  : the family index of the arrays to redistribute
c int DISTTYPE_BRA : the distribution type (cf. disttype.h) of the bra (rows)
c int DISTTYPE_KET : the distribution type (cf. disttype.h) of the ket (cols)
c logical BUPDATE : a behavior flag
c                   = T;        update /sympop/isytyp(*,IFAMNDX)
c                   = F; do not update ...
c                   Note: This flag is so important that it is ignored
c                         and assumed to be true.





c#define _DEBUG_ACES_FAMILY_REDIST

      subroutine aces_family_redist(irp_diag,iFamNdx,
     &                              DistType_bra,DistType_ket,
     &                              bUpdate)
      implicit none

c ARGUMENTS
      integer irp_diag, iFamNdx, DistType_bra, DistType_ket
      logical bUpdate

c INTERNAL VARIABLES
      integer irp_bra, irp_ket
      integer nRows, nCols, nWords
      integer iNumNew, iNumMax
      integer iRec, iWrdNdx
      integer quiktmp(10)
      integer iTmp

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
c filspc.com : begin

c This common block contains the dimensions of the physical records used by the
c MOIO storage files.

c iprcln : the    byte-length of a physical record
c iprcwd : the integer-length of a physical record
c iprcrl : the    recl-length of a physical record

      integer         iprcln, iprcwd, iprcrl
      common /filspc/ iprcln, iprcwd, iprcrl
      save   /filspc/

c filspc.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
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


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end




c ----------------------------------------------------------------------

      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert irp_diag is in [1,nirrep]
      if ((irp_diag.lt.1).or.(nirrep.lt.irp_diag)) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   irp_diag = ',irp_diag
         print *, '   nirrep   = ',nirrep
         iTmp = 1
      end if
c   o assert iFamNdx is properly bound
      if ((iFamNdx.lt.1).or.(500.lt.iFamNdx)) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   iFamNdx = ',iFamNdx
         print *, '   max ndx = ',500
         iTmp = 1
      end if
c   o assert distribution types are in [1,22]
      if ((DistType_bra.lt.1).or.(22.lt.DistType_bra).or.
     &    (DistType_ket.lt.1).or.(22.lt.DistType_ket)    ) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   DistType_bra = ',DistType_bra
         print *, '   DistType_ket = ',DistType_ket
         print *, '   max type     = ',22
         iTmp = 1
      end if
c   o assert bUpdate is .true. (even though it's ignored)
      if (.not.bUpdate) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   bUpdate = ',bUpdate
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)

c ----------------------------------------------------------------------

c   o process the first row/group/irrep (peeled off the do loop)
      irp_bra = dirprd(1,irp_diag)
      nRows   = irpdpd(irp_bra,DistType_bra)
      nCols   = irpdpd(1,      DistType_ket)
      iNumNew = nRows*nCols
      iNumMax = moiomxsz(1,iFamNdx)
      moiosz(1,iFamNdx) = nRows
      moiods(1,iFamNdx) = nCols

c   o process 'the others'
      if (nirrep.ne.1) then
         nWords = nRows*nCols*iintfp
         quiktmp(1) = nWords
         iWrdNdx = moiowd(1,iFamNdx) + nWords
         iTmp    = (iWrdNdx-1)/iprcwd
         iRec    = moio(1,iFamNdx) + iTmp
         iWrdNdx = iWrdNdx         - iTmp*iprcwd
         do irp_ket = 2, nirrep
            irp_bra = dirprd(irp_ket,irp_diag)
            nRows   = irpdpd(irp_bra,DistType_bra)
            nCols   = irpdpd(irp_ket,DistType_ket)
            iTmp    = nRows*nCols
            iNumNew = iNumNew + iTmp
            nWords  = iTmp*iintfp
            iNumMax = iNumMax + moiomxsz(irp_ket,iFamNdx)
            moiosz(irp_ket,iFamNdx) = nRows
            moiods(irp_ket,iFamNdx) = nCols
            moio  (irp_ket,iFamNdx) = iRec
            moiowd(irp_ket,iFamNdx) = iWrdNdx
            quiktmp(irp_ket) = nWords
            iWrdNdx = iWrdNdx + nWords
            iTmp    = (iWrdNdx-1)/iprcwd
            iRec    = iRec    + iTmp
            iWrdNdx = iWrdNdx - iTmp*iprcwd
         end do
         if (quikget(1,iFamNdx).ne.0) then
            do irp_ket = 2, nirrep
               quikget(irp_ket,iFamNdx) =   quikget(irp_ket-1,iFamNdx)
     &                                    + quiktmp(irp_ket-1)
            end do
         end if
c     end if (nirrep.ne.1)
      end if

cYAU - I don't think we should modify moiomxsz. Those elements are set in
c      aces_list_touch which also moves the free-space pointers.

c   o make sure programmers are not creating new data
      if (iNumNew.gt.iNumMax) then
      end if

c   o update the type lookup table
      isytyp(1,iFamNdx) = DistType_bra
      isytyp(2,iFamNdx) = DistType_ket

c   o mark MOIO as modified
      bIOMod = .true.

      return
c     end subroutine aces_family_redist
      end

