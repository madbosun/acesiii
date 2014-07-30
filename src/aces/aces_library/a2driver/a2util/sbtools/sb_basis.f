
c This routine reads the MOL file and initializes data in /MOL/.
c It should be called only by sb_com_mol(), but its main purpose is
c to recast icore() addresses into convenient variables.
c
c Determines totprim, maxshlprim, maxshlorb, and maxshell.  It also
c fills the arrays nprimatom, nshellatom, and naoatom.

      subroutine sb_basis(comppopv,compmemb,nangatom,nshellatom,
     &                    nprimatom,naoatom)
      implicit none



c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)




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



      integer comppopv(*),compmemb(*),nangatom(*),nshellatom(*),
     &        nprimatom(*),naoatom(*)

      integer
     &    iorbit,nshell,ishell,nsubshell(maxangshell),isub,natomprim,
     &    natomao,natomang,ishellang,ishellprim,ishellao,isubprim,
     &    isubao,nline,icount,ieqatm,j,iatom,iunit
      character*80 line,molfil

      integer maxorbit
      parameter (maxorbit=1000)
      integer orbitpop(maxorbit)

      call callstack_push('SB_BASIS')

c See mol.com for an explanation of this.  It will disappear soon.
      if (memknown.ne.1) then
         call getrec(0,'JOBARC','COMPPOPV',j,orbitpop)
         if (j.gt.maxorbit) then
            print *, '@SB_BASIS: hard limit exceeded; increase maxorbit'
            print *, '           dim(COMPPOPV) = ',j
            call errex
         end if
         call getrec(1,'JOBARC','COMPPOPV',natomsx,orbitpop)
      end if

      iunit=3
      call gfname('MOL',molfil,j)
      open(unit=iunit,file=molfil(1:j),form='formatted',status='old')
      rewind(iunit)
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

      icount=0
      totprim=0
      maxshlprim=0
      maxshlao=0
      maxshlorb=0
      maxshell=0
      maxangmom=0
      do iorbit=1,compnorb
         read(iunit,'(a25,12i5)') line,nshell,
     &                            (nsubshell(ishell),ishell=1,nshell)
         read(iunit,'(a)') line
         natomprim=0
         natomao=0
         natomang=0
         do ishell=1,nshell
c           ishellang= 1,3,6,10,...   (1s, 3p, 6d, 10f, ...)
            ishellang=rshift(ishell*(ishell+1),1)
            ishellprim=0
            ishellao=0
            do isub=1,nsubshell(ishell)
               read(iunit,'(2i5)') isubprim,isubao
               ishellprim=ishellprim+isubprim
               ishellao=ishellao+isubao

               nline=(isubao-3)/4
               if ((isubao-3).gt.(nline*4)) nline=nline+1
               nline=(nline+1)*isubprim
               do j=1,nline
                  read(iunit,'(a)') line
               end do
            end do
            natomao=natomao+ishellang*ishellao
            natomprim=natomprim+ishellang*ishellprim
            natomang=natomang+ishellang
            if (ishellprim.gt.maxshlprim) maxshlprim=ishellprim
            if (ishellao.gt.maxshlao) maxshlao=ishellao
            if ((ishellprim*ishellao).gt.maxshlorb)
     &         maxshlorb=ishellprim*ishellao
c        end do ishell=1,nshell
         end do
         if (natomang.gt.maxshell) maxshell=natomang
         if (nshell.gt.maxangmom) maxangmom=nshell
         if (memknown.ne.0) then
            do ieqatm=1,comppopv(iorbit)
               icount=icount+1
               iatom=compmemb(icount)
               nprimatom(iatom)=natomprim
               nshellatom(iatom)=natomang
               naoatom(iatom)=natomao
               nangatom(iatom)=nshell
            end do
         end if
         totprim=totprim+natomprim*orbitpop(iorbit)
c     end do iorbit=1,compnorb
      end do
      
      close(iunit)
      call callstack_pop
      return
      end

