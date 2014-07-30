
c This routine gets all of the standard information from the JOBARC file
c for the molecule.

      subroutine sb_com_mol
      implicit none















c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler


























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






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





c Two different "stacks" of memory will be available in the course of an
c Aces3 calculation.  Integer memory will be stored in icore.  Real memory
c will be stored in dcore.
c
c Memory allocation is done in one of two ways, dynamic or nondynamic.  This
c is controlled by the parameter dynmem.  If dynmem is 1, dynamic memory
c is used.  Otherwise, nondynamic memory is used.
c
c Two additional parameters nondynimem and nondyndmem control how much
c nondynamic memory is allocated initially in the icore and dcore stacks.
c If dynmem is 1, both of these parameters SHOULD be 1.

c Dynamic memory allocation is done in a fairly straighforward way.  The
c program runs through everything twice.  The first time is just to determine
c how much dynamic memory is required.  The second time is use to actually
c allocate the memory and use it.  In some cases, it may be difficult to
c determine how much memory to use in advance in each of the stacks.  To
c aid this, there is also the option of only running through a program a
c single time.  When this happens, both stacks are allocated initially and
c used throughout the program.  If the memory requirement for either stack
c is exceeded, the program crashes.
c
c Two parameters which are used in either case are memknown and maxicore.
c Both MUST be set in the calling program.
c
c   memknown : 0 if a dummy run is being done to determine memory usage
c              1 if the memory usage is known from a previous dummy run
c             -1 if no dummy run is done (only one run through the program)
c   maxicore : the maximum amount of icore to allocate (this MAY be
c              adjusted IF a dummy run detects that more is needed, but in
c              the case of a single run, this will be rigidly applied)
c
c This include file is required whenever calls to any of the memory functions
c (setptr and relptr) are called or when part of the execution depends on
c whether the first or second run is being performed.
c
c ***NOTE***
c It is important that icore be aligned on a floating point boundary.  One
c way which seems to insure this is to have icore the FIRST element in the
c common block.  So, make sure that no integer is ever inserted before icore
c in the common block declaration.

      integer dynmem,nondynimem,nondyndmem
      parameter (dynmem=1)
      parameter (nondynimem=1)
      parameter (nondyndmem=1)

      double precision   dcore(nondyndmem)
      common /dcore_com/ dcore
      save   /dcore_com/

      integer    kscore(nondynimem)
      common /kscore_mem / kscore

      integer            memknown,maxicore
      common /kscore_com/ memknown,maxicore
      save   /kscore_com/





      integer maxirrep,num2comb,max2comb
      parameter (maxirrep=8)
      parameter (num2comb=22)
      parameter (max2comb=25)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer        nirrep, numbasir(8),
     &               irpsz1(36),irpsz2(28),irpds1(36),irpds2(56),
     &               old_irpoff(9), irrorboff(9), dirprd(8,8),
     &               old_iwoff1(37), old_iwoff2(29),
     &               inewvc(maxbasfn), idxvec(maxbasfn),
     &               irrtrilen(9), irrtrioff(8),
     &               irrsqrlen(9), irrsqroff(8)
      common /symm2/ nirrep, numbasir,
     &               irpsz1,    irpsz2,    irpds1,    irpds2,
     &               old_irpoff,    irrorboff,    dirprd,
     &               old_iwoff1,     old_iwoff2,
     &               inewvc,           idxvec,
     &               irrtrilen,    irrtrioff,
     &               irrsqrlen,    irrsqroff
      save   /symm2/

      integer             occup(8,2),totocc(2),totocca,totoccb,
     &                    maxirrtri,maxirrsqr,irrtritot,irrsqrtot
      common /sym_ks_com/ occup,     totocc,   totocca,totoccb,
     &                    maxirrtri,maxirrsqr,irrtritot,irrsqrtot
      save   /sym_ks_com/




c This common block contains the molecule and basis set information
c read in from the JOBARC file.  Since much of this information is
c used in a large number of modules, and since most of the information
c is relatively small compared to the other things held in memory,
c a large percentage of the data stored in the JOBARC file is stored
c here, even though some modules will not use all of it.

c   maxangshell - The maximum number of angular momentum shells.  Since this
c                 is used VERY infrequently, set it high enough to never
c                 cause a problem.
c   spinc(2)    - The characters 'A' and 'B' (useful for alpha/beta labels)
c   natoms      - The number of atoms in the Z-matrix (including X/GH).  After
c                 remove is called, natoms becomes equivalent to nrealatm.
c   natomsx     - The number of atoms in the Z-matrix (including X/GH).  This
c                 does not change.
c   nrealatm    - The number of atoms in the Z-matrix (including GH).
c   naobasfn    - The number of AOs in the basis
c   nbastot     - The number of symmetry adapted orbitals in the basis (the AO
c                 basis may be larger than the SO basis if spherical orbitals
c                 are used since harmonic contaminants are deleted)
c   linear      - 1 if the molecule is linear
c   orientmt    - 3x3 matrix which relates the computational and canonical
c                 orientations
c   nucrep      - Nuclear repulsion energy in a.u.
c   nmproton    - Number of protons in the molecule.
c
c   compptgp    - Point group
c   fullptgp    -
c   compordr    - Order of the point group
c   fullordr    -
c   compnirr    - Number of irreps in the point group
c   fullnirr    -
c   compnorb    - Number of unique atoms (orbits) in the point group
c   fullnorb    -
c   c1symmet    - 1 if the molecule is C1 symmetry
c   nirrep      - The same as compnirr (since nirrep is used so commonly,
c                 this is included for conveniance)
c                 ***NOTE*** nirrep is read in twice and is stored in /sym/
c                            so it is not actually included here
c   totprim     - Total number of primitive functions in the molecule
c   maxshlprim  - Largest number of primitives in a single shell
c   maxshlao    - Largest number of AOs in a single shell
c   maxshlorb   - Largest number of primitive orbitals (primitive functions
c                 times the number of AOs) in a single shell
c   maxangmom   - Largest angular momentum for any atom
c   maxshell    - Larges number of angular momentum shells for any atom
c   noccorb(2)  - The number of alpha and beta occupied orbitals
c   nvrtorb(2)  - The number of alpha and beta virtual orbitals

c The parameter maxorbit is needed because of how dynamic memory is used.
c Two runs of the program are needed.  The first to calculate memory usage,
c the second to use it.  In order to calculate totprim, we have to know the
c orbit population vector (the number of each type of atom).  BUT, this is
c stored in dynamic memory since we do not know how long this vector is.
c In the future, joda or vmol will write this information to JOBARC, and
c this problem will disappear.  In the meantime, we have to introduce a
c genuine limit on the size of the molecule.  It may have no more than
c maxorbit sets of unique atoms.  This limit is ONLY used in the subroutine
c basis, so it probably will disappear when the information in the MOL file
c is put in JOBARC.
c    maxorbit   - the number of symmetry unique atoms

c The following are pointers to real arrays
c
c   zatommass(natoms)  - Atomic mass of all atoms (X=0.0, GH=100.0)
c   zcoord(3,natoms)   - Coordinates of all atoms (computational orientation)
c   zalpha(totprim)    - The alpha for each primitive function
c   zprimcoef(totprim,naobasfn)
c                      - The primitive to AO coefficients
c
c The following are pointers to integer arrays
c
c   patomchrg(natoms)  - Atomic number of all atoms (X=0, GH=110)
c   pfullclss(fullordr)- Class type vector
c   pcompclss(compordr)-
c   pfullpopv(natoms)  - Number of atoms in each orbit
c   pcomppopv(natoms)  -
c   pfullmemb(natoms)  - Atoms sorted by point group orbits
c   pcompmemb(natoms)  -
c   pnprimatom(natoms) - Number of primitive functions for each atom
c   pnshellatom(natoms)- Number of different angular momentum shells for each
c                        atom (takes on values of 1,4,9,16, etc.)
c   pnangatom(natoms)  - The number of different angular momentum for each
c                        atom (takes on values of 1,2,3,4, etc.)
c   pnaoatom(natoms)   - Number of AOs for each atom
c   pnshellprim(maxshell,natoms)
c                      - The number of primitive functions in each shell
c                        of each atom
c   pnshellao(maxshell,natoms)
c                      - The number of AOs in each shell of each atom
c   pprimoff(maxshell,natoms)
c   paooff(maxshell,natoms)
c                      - The primcoef matrix is a block diagonal matrix.
c                        Each shell of each atom has a block.  If you have
c                        a list of all primitive functions, pprimoff(ishell,
c                        iatom) tells the location of the first primitive
c                        function in the block (ishell,iatom) and paooff
c                        contains similar information for the AOs.
c
c ***NOTE***  Because joda stores pfullpopv/pcomppopv as size natoms, we
c             do to, but they should be of size fullnorb/compnorb.  The
c             first ones have real values.  The remaining ones are 0.

      double precision orientmt(3,3),nucrep
      integer natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb(2),
     &    nvrtorb(2),maxshlao,maxangmom,natomsx
      integer patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,
     &    zalpha,zprimcoef,pprimoff,paooff,pnangatom

      common /mol_com/ orientmt,nucrep,
     &    natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb,
     &    nvrtorb,maxshlao,maxangmom,natomsx,
     &    patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,zalpha,
     &    zprimcoef,pprimoff,paooff,pnangatom
      save   /mol_com/

      character*4 compptgp,fullptgp
      character*1 spinc(2)
      common /molc_com/ compptgp,fullptgp,spinc
      save   /molc_com/



      integer setptr,i
      character*8 szInt8
      integer    zccoeff,pprims,ppriml
      call callstack_push('SB_COM_MOL')

      call getrec(1,'JOBARC','NATOMS',1,natoms)
      natomsx=natoms
      if (memknown.ne.1) then
         call getrec( 1,'JOBARC','NREALATM',1,         nrealatm)
         call getrec( 1,'JOBARC','NUCREP  ',iintfp,    nucrep)
         call getrec(-1,'JOBARC','LINEAR  ',1,         linear)
         call getrec(-1,'JOBARC','ORIENTMT',3*3*iintfp,orientmt)
         call getrec(-1,'JOBARC','NMPROTON',1,         nmproton)
         call getrec(-1,'JOBARC','NOCCORB ',2,         noccorb)
         call getrec(-1,'JOBARC','NVRTORB ',2,         nvrtorb)
      end if

      patomchrg=setptr(1,0,natoms)
      zatommass=setptr(1,1,   natoms)
 
      zcoord   =setptr(1,1,   3*natoms)

      if (memknown.ne.0) then
         call getrec(1,'JOBARC','ATOMCHRG',natoms,kscore(patomchrg))
         call getrec(1,'JOBARC','ATOMMASS',natoms*iintfp,
     &               dcore(zatommass))
         call getrec(1,'JOBARC','COORD',3*iintfp*natoms,dcore(zcoord))
      end if

      if (memknown.ne.1) then
c        Reading strings!!!
         call getrec(1,'JOBARC','COMPPTGP',iintfp,szInt8)
         compptgp(1:4)=szInt8(1:4)
         call getrec(1,'JOBARC','FULLPTGP',iintfp,szInt8)
         fullptgp(1:4)=szInt8(1:4)
         c1symmet=0
         if (compptgp.eq.'C1  ') c1symmet=1

         compnirr=nirrep
         call getrec( 1,'JOBARC','FULLNIRR',1,fullnirr)
         call getrec( 1,'JOBARC','COMPNORB',1,compnorb)
         call getrec( 1,'JOBARC','FULLNORB',1,fullnorb)
         call getrec(-1,'JOBARC','COMPORDR',1,compordr)
         call getrec(-1,'JOBARC','FULLORDR',1,fullordr)

         spinc(1)='A'
         spinc(2)='B'
      end if

      pfullclss=setptr(1,0,fullordr)
      pcompclss=setptr(1,0,compordr)
      pfullpopv=setptr(1,0,natoms)
      pcomppopv=setptr(1,0,natoms)
      pfullmemb=setptr(1,0,natoms)
      pcompmemb=setptr(1,0,natoms)
      if (memknown.ne.0) then
         call getrec(-1,'JOBARC','FULLCLSS',fullordr,kscore(pfullclss))
         call getrec(-1,'JOBARC','COMPCLSS',compordr,kscore(pcompclss))
         call getrec( 1,'JOBARC','FULLPOPV',natoms,  kscore(pfullpopv))
         call getrec( 1,'JOBARC','COMPPOPV',natoms,  kscore(pcomppopv))
         call getrec( 1,'JOBARC','FULLMEMB',natoms,  kscore(pfullmemb))
         call getrec( 1,'JOBARC','COMPMEMB',natoms,  kscore(pcompmemb))
      end if

c   o remove dummy atoms
      if (memknown.ne.0) then
         call sb_removedummy(dcore(zatommass),dcore(zcoord),
     &                       kscore(patomchrg),kscore(pfullmemb),
     &                       kscore(pcompmemb),natomsx)
      end if
      natoms=nrealatm

c ----------------------------------------------------------------------
c sb_getbas()

      if (memknown.ne.1) then
         call getrec(1,'JOBARC','NAOBASFN',1,naobasfn)
         call getrec(1,'JOBARC','NBASTOT ',1,nbastot)
      end if

      pnprimatom =setptr(1,0,natoms)
      pnshellatom=setptr(1,0,natoms)
      pnangatom  =setptr(1,0,natoms)
      pnaoatom   =setptr(1,0,natoms)


      call sb_basis(kscore(pcomppopv),kscore(pcompmemb),
     &              kscore(pnangatom),kscore(pnshellatom),
     &              kscore(pnprimatom),kscore(pnaoatom))

      pnshellprim=setptr(1,0,maxshell*natoms)
      pnshellao  =setptr(1,0,maxshell*natoms)
      pprimoff   =setptr(1,0,maxshell*natoms)
      paooff     =setptr(1,0,maxshell*natoms)

      zalpha     =setptr(1,1,totprim)
      zprimcoef  =setptr(1,1,totprim*naobasfn)

      zccoeff   =setptr(1,1,totprim*naobasfn)
      pprims   =setptr(1,0,naobasfn)
      ppriml   =setptr(1,0,naobasfn)

      if (memknown.ne.0) then
         call sb_prim(kscore(pcomppopv),kscore(pcompmemb),
     &                kscore(pnshellatom),
     &                kscore(pnshellprim),kscore(pnshellao),
     &                kscore(pprimoff),kscore(paooff),
     &                dcore(zalpha),dcore(zprimcoef),
     &     dcore(zccoeff),kscore(pprims),kscore(ppriml))
      end if

c ----------------------------------------------------------------------

      call callstack_pop
      return
      end

