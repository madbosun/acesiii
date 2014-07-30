C
cjp
cjp in bwcc it computes diagonal elements of the effective hamiltonian,
cjp not the energy, which is obtained from diagonalization of that Heff
cjp after also the offdiagonal elements are computed
cjp
C  DRIVER FOR THE  CALCULATION OF THE CORRELATION ENERGY FOR A GIVEN SET
C  OF AMPLITUDES
C
C  ARGUMENTS :  ICORE ..... ICORE ARRAY
C               MAXCOR .... DIMENSION OF ICORE
C               NLIST2 .... OFFSET OF T2 LIST ON MOINTS (WITH RESPECT TO
C                            TYPE)
C               NLIST1 .... OFFSET OF T1 LISTS ON MOINTS (WITH RESPECT TO
C                            SPIN TYPE)
C               ECORR ..... RETURNS THE CORRELATION ENERGY FOR ALL SPIN CASES
C               ETOT .....  RETURNS TOTAL CORRELATION ENERGY
C               ETOTT2 ...  RETURNS LINEAR CONTRIBUTION TO
C                            THE CORRELATION ENERGY
C               IUHF .....  IUHF FLAG
C               IPRINT ...  PRINT FLAG
C
      SUBROUTINE CMPENG(ICORE,MAXCOR,NLIST2,NLIST1,ECORR,ETOT,ETOTT2,
     &                  IUHF,IPRINT)
      IMPLICIT INTEGER (A-Z)
      LOGICAL TAU,NONHF
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1
      CHARACTER*2 SPCASE(3)
      DOUBLE PRECISION E,ETOT,FACTOR,ECORR(3),ESPIN,ET2,ETOTT2,
     &                 ESING,SDOT
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD
      COMMON /NHFREF/NONHF
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DIMENSION I0T(2),I0F(2)
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA SPCASE /'AA','BB','AB'/
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



C
c Nevin insured that mxcor is even for alignment
c      MXCOR=MAXCOR
      MXCOR=MAXCOR - MOD(MAXCOR,2)

cYAU: The one-particle entities were moved from the end of icore to the front.
      IFREE=1

      TAU=.FALSE.
      IF((METHOD.GT.9.AND.SING1.AND.METHOD.NE.21.AND.METHOD.NE.23)
     &    .OR.NONHF.AND.SING1)THEN
C
C   ALLOCATE MEMORY FOR T1 AMPLITUDES
C
cYAU       I0T(1)=MXCOR+1-NT(1)*IINTFP
cYAU       MXCOR=MXCOR-NT(1)*IINTFP
       I0T(1)=IFREE
       IFREE=IFREE+NT(1)*IINTFP
       CALL GETLST(ICORE(I0T(1)),1,1,1,1+NLIST1,90)
       IF(IUHF.EQ.0) THEN
        I0T(2)=I0T(1)
       ELSE
cYAU        I0T(2)=I0T(1)-NT(2)*IINTFP
cYAU        MXCOR=MXCOR-NT(2)*IINTFP
        I0T(2)=IFREE
        IFREE=IFREE+NT(2)*IINTFP
        CALL GETLST(ICORE(I0T(2)),1,1,1,2+NLIST1,90)
       ENDIF
       TAU=.TRUE.
C
C   FOR NON HF REFERENCES ALLOCATE MEMORY FOR f(a,I)
C
       IF(NONHF) THEN
cYAU        I0F(1)=I0T(2)-NT(1)*IINTFP
cYAU        MXCOR=MXCOR-NT(1)*IINTFP
        I0F(1)=IFREE
        IFREE=IFREE+NT(1)*IINTFP
        CALL GETLST(ICORE(I0F(1)),1,1,1,3,93)
        IF(IUHF.EQ.0) THEN
         I0F(2)=I0F(1)
        ELSE
cYAU         I0F(2)=I0F(1)-NT(2)*IINTFP
cYAU         MXCOR=MXCOR-NT(2)*IINTFP
         I0F(2)=IFREE
         IFREE=IFREE+NT(2)*IINTFP
         CALL GETLST(ICORE(I0F(2)),1,1,1,4,93)
        ENDIF
       ENDIF
      ENDIF
C
      ETOT=0.D0
      ETOTT2=0.D0
      FACTOR=1.D0
      IF(IUHF.EQ.0)FACTOR=2.D0
      DO 10 ISPIN=1,IUHF+1
       LISTT=43+ISPIN
       ESPIN=0.D0
C
C  FOR NON HF REFERENCE FUNCTION ADD HERE THE T(I,A) f(I,A) CONTRIBUTION
C  THERE ARE AA AND BB CONTRIBUTIONS
C
       IF(NONHF.AND.SING1) THEN
        ESING=SDOT(NT(ISPIN),ICORE(I0T(ISPIN)),1,ICORE(I0F(ISPIN)),1)
        ESPIN=ESPIN+ESING
        ETOT=ETOT+FACTOR*ESING
        ETOTT2=ETOTT2+FACTOR*ESING
       ENDIF
C
C  THE TAU(IJ,AB) <IJ//AB> CONTRIBUTION TO THE ENERGY (AA AND BB SPIN
C  CASES
C
       DO 100 IRREP=1,NIRREP
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
cYAU        I001=1
        I001=IFREE
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(I004.LT.MXCOR) THEN
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &             ICORE(I002),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ISPIN,TAU,IRREP,POP(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2=ETOTT2+FACTOR*ET2
        ETOT=ETOT+FACTOR*E
        ESPIN=ESPIN+E
        ENDIF
100    CONTINUE
        ECORR(ISPIN)=ESPIN
       IF(IPRINT.NE.0)WRITE(*,80)SPCASE(ISPIN),ESPIN
80     FORMAT(T3,' The ',A2,' contribution to the correlation ',
     &        'energy is: ',F12.7,' a.u.')
10    CONTINUE
C
C  THE TAU(Ij,Ab) <Ij//Ab> CONTRIBUTION TO THE ENERGY (SPIN CASE AB)
C
      ESPIN=0.D0
      DO 200 IRREP=1,NIRREP
       DISSYT=IRPDPD(IRREP,ISYTYP(1,46))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,46))
       IF(MIN(NUMSYT,DISSYT).NE.0) THEN
cYAU        I001=1
        I001=IFREE
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(I004.LT.MXCOR) THEN
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &              ICORE(I002),ICORE(I0T(1)),ICORE(I0T(2)),3,TAU,
     &              IRREP,POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),
     &              ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2=ETOTT2+ET2
        ETOT=ETOT+E
        ESPIN=ESPIN+E
       ENDIF
200   CONTINUE
      ECORR(3)=ESPIN
      IF(IPRINT.NE.0)WRITE(*,80)SPCASE(3),ESPIN
      IF(IPRINT.NE.0)WRITE(*,81)ETOT
81    FORMAT(T3,' The total correlation energy is ',F15.12,' a.u.')
cjp
cjp add vacuum energy of (iref) and store the result in diagonal element
cjp of the effective hamiltonian
cjp
      if (isbwcc) heff(iref,iref)=etot+hdiagcontr(iref)
      RETURN
      END
