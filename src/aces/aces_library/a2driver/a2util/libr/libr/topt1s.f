      SUBROUTINE TOPT1S(ISPIN,ILIST,NTOP,NT1SIZ,NOCC,
     &                  NVRT,NSMSZ,T,TOPT1,
     &                  SYVEC,ITOPT1,I,A,TYPE)
C
C THIS ROUTINE PICKS UP A T1 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT1 AND THEIR ASSOCIATED INDICES
C  IN THE I AND A VECTORS.
C
C INPUT:
C       ISPIN - THE LIST SUBTYPE
C       ILIST - THE LIST NUMBER
C       NTOP  - THE NUMBER OF LARGEST AMPLITUDES (BY ABSOLUTE VALUE)
C               WHICH ARE TO BE FOUND
C       NT1SIZ- THE TOTAL SIZE OF THE T1 VECTOR
C       NOCC  - THE NUMBER OF OCCUPIED ORBITALS OF THIS SPIN TYPE
C       NVRT  - THE NUMBER OF VIRTUAL ORBITALS OF THIS SPIN TYPE
C       NSMSZ - THE TOTAL SIZE OF THE SYMMETRY VECTOR (NOCC*NVRT)
C       SYVEC - THE SYMMETRY VECTOR FOR THE A,I OR a,i LIST
C       TYPE  - ???????????
C
C OUTPUT:
C       TOPT1 - THE NTOP LARGEST ELEMENTS IN THE T1 VECTOR
C               SORTED BY ABSOLUTE VALUE
C       I     - THE I INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       A     - THE A INDICES CORRESPONDING TO THE VALUES IN TOPT1
C
C SCRATCH:
C
C       T     - USED TO HOLD THE SYMMETRY-PACKED T VECTOR
C       ITOPT1- USED TO HOLD THE OFFSETS CORRESPONDING TO THE
C               ELEMENTS IN TOPT1
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A(NTOP),SYVEC(NSMSZ)
      CHARACTER*2 SPCASE(3)
      CHARACTER*1 TYPE
      DIMENSION T(NT1SIZ),TOPT1(NTOP),ITOPT1(NTOP),I(NTOP)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
cjp





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



      DATA SPCASE /'AA','BB','AB'/
C
C THIS ROUTINE PICKS UP A T1 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT1 AND THEIR ASSOCIATED INDICES
C  IN THE I AND A VECTORS.
C
C
C STATEMENT FUNCTIONS
C
      IGETI(IX,NUMA)=1+(IX-1)/NUMA
      IGETA(IX,NUMA)=IX-(IGETI(IX,NUMA)-1)*NUMA
      NNM1O2(IX)=(IX*(IX-1))/2
C
      IF(NT1SIZ.EQ.0)RETURN
C
      CALL GETLST(T,1,1,1,ISPIN,ILIST)
C
C COMPUTE NORM OF T1 VECTOR
C
      X=SNRM2(NT1SIZ,T,1)
      T1NORM=X
C
C PUT NTOP LARGEST ELEMENTS INTO TOPT1 AND THEIR ASSOCIATED INDICES INTO
C   ITOPT1.
C
      CALL SCANVC(T,TOPT1,ITOPT1,NTOP,NT1SIZ)
      DO 10 IRANK=1,NTOP
       IUNPK=SYVEC(ITOPT1(IRANK))
       A(IRANK)=IGETA(IUNPK,NVRT)+NOCCO(ISPIN)
       I(IRANK)=IGETI(IUNPK,NVRT)
10    CONTINUE
cjp
      if(isbwcc) then
        write(*,201)type,SPCASE(ISPIN),iref
        else
      WRITE(*,200)TYPE,SPCASE(ISPIN)
        endif
200   FORMAT(T3,' Largest ',A1,'1 amplitudes for spin case ',A2,':')
201   format(T3,' Largest ',A1,'1 amplitudes for spin case ',A2,
     +   ' of reference no.: ',i3)
      IF(ISPIN.EQ.3)WRITE(*,299)
      WRITE(*,300)
299   FORMAT(3X,3(' ',3X,'_',3X,' ',3X,'_',13X))
300   FORMAT(3X,3('i',3X,' ',3X,'a',3X,' ',13X))
      WRITE(*,400)
400   FORMAT(77('-'))
      if(isbwcc) then
      if(bwgossip) then
        write(*,501) (irank,i(irank),a(irank),
     &         invhnum(i(irank),iref,ispin),
     &        invpnum(a(irank)-nocco(ispin),iref,ispin),
     &         topt1(irank),irank=1,ntop)
        else
         write(*,500)(invhnum(i(irank),iref,ispin),
     &      invpnum(a(irank)-nocco(ispin),iref,ispin),
     &              topt1(irank),irank=1,ntop)
        endif
      else
      WRITE(*,500)(I(IRANK),A(IRANK),
     &                 TOPT1(IRANK),IRANK=1,NTOP)
      endif
      WRITE(*,400)
      WRITE(*,600)TYPE,SPCASE(ISPIN),NT1SIZ,T1NORM
      WRITE(*,400)
600   FORMAT(T3,' Norm of ',A1,'1',A2,' vector (',I9,
     &        ' symmetry allowed elements):',F14.10,'.')
500   FORMAT((2('[',I3,1X,3X,1X,I3,1X,3X,']',F8.5,1X),'[',I3,1X,
     &3X,1X,I3,1X,3X,']',F8.5))
cjp
501    format(i4,':',i4,1x,i4,':true indices:',i4,1x,i4,' : ',f12.8)
      RETURN
      END
