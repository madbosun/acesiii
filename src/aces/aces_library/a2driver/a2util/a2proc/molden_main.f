
c This program writes a molden format file containing the coordinates,
c basis functions for each center, molecular orbitals, frequencies,
c frequency coordinates, and normal mode vectors.
c
c Basis functions are limited to gaussian which is a constraint from aces2.
c
c the following files need to be present:
c	JOBARC
c	JAINDX
c	ZMAT
c	GENBAS
c       for frequencies NORMCO
c 
c Ken Wilson March 1998































































































































































































































































































































































































































































































































      subroutine molden_main (args,dimargs)

      implicit double precision (a-h,o-z)
C
C Watson added
C 
c ARGUMENT LIST
      integer dimargs
      character*80 args(dimargs)
      
      logical LNAT_ORBS, dens_diff
      character*8 DLABEL
      character*2 iroot

      parameter (mxangmom=7)
      parameter (mxcoef=30)


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end







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



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
      character*32 szFile
      logical bExist
      dimension iangmom(mxangmom)
c
      iuhf = 1
      if (iflags(11).eq.0) iuhf = 0
      maxcor=icrsiz
      iunit =10
      if (iunit.ne.6) then
         szFile = 'MOLDEN.INPUT'
         inquire(file=szFile,exist=bExist)
         if (bExist) call f_remove(szFile)
         open(unit=iunit,file=szFile,form='formatted',err=2)
 1       goto 3
 2       print *, '@MOLDEN_MAIN: could not create/open ',szFile
         call errex
 3       continue
      end if
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      iexx=0
c      iexx=1
c      call nl_init('tomolden',ierr,.true.)
c      if(ierr.eq.1)then
c         write(*,*)' *tomolden namelist not found in ZMAT file'
c         call errex
c      end if
c      call nl_str('EXX',0,iexx)
c      call nl_term
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c geometry
c
      write(iunit,'(a)')'[Molden Format]'
      write(iunit,'(a)')'[Atoms] Angs'

      call getrec(20,'JOBARC','NATOMS',1,natoms)

      iatchrg=i0
      icoord=iatchrg+natoms+mod(natoms,2)
      ifuchrg=icoord+3*natoms*iintfp
      ifucoord=ifuchrg+natoms+mod(natoms,2)
      inext=ifucoord+3*natoms*iintfp

      if(inext-i0.gt.maxcor)call insmem('rdgeom-f',inext-i0,maxcor)

      call A2rd_geom(natoms,icore(iatchrg),icore(icoord),
     &               icore(ifuchrg),icore(ifucoord),iunit,.true.)
c--------------------------------------------------------------------------
c basis set
c
      call getrec(0,'JOBARC','SCFENEG ',Ilength,Ijunk)
      If (Ilength .LE. 0) Then
         if (iunit.ne.6) then
             write(iunit,*)
             close(iunit)
             print *, '@MOLDEN_MAIN: successfully created ',szFile
         end if
         Return
      Endif

      write(iunit,'(a)')'[GTO]'

      iexp=inext 
      icoef=iexp+mxangmom*mxcoef*iintfp 
      inext=icoef+mxangmom*mxcoef**2*iintfp 

      if(inext-i0.gt.maxcor)call insmem('rdbasis-f',inext-i0,maxcor)

      call A2rd_basis(natoms,icore(iatchrg),icore(iexp),
     &                icore(icoef),nshells,iangmom,iunit,.true.)
c--------------------------------------------------------------------------
c orbitals from exact calculations
c
      if((iflags(54).ne.3).and.
     &   (iflags2(138).ne.3)) then
   
        call getrec(0,'JOBARC','SCFENEG ',Ilength,Dtmp)
        if (ilength .gt. 0) then
           write(iunit,*)
           write(iunit,'(a)')'[MO]'
        else
           if (iunit.ne.6) then
              close(iunit)
              print *, '@MOLDEN_MAIN: successfully created ',szFile
           end if
            Return
        Endif
C
        call getrec(-1,'JOBARC','NAOBASFN',1,nao)
        call getrec(-1,'JOBARC','NBASTOT',1,nmo)
C
C Watson, added
C
        iroot = "00"
        dens_diff=.FALSE.
        LNAT_ORBS = .FALSE.
        iLst = linblnk(args(1))
        IF (args(1)(1:iLst) .EQ. "nlo") LNAT_ORBS = .TRUE.

        iLst = linblnk(args(2))
        IF (iLst .NE. 0) THEN
          DLABEL = args(2)(1:iLst)
          IF (DLABEL(1:iLst) .EQ. "ediff") THEN
            iLst = linblnk(args(3))
            dens_diff=.TRUE.
            IF (iLst .EQ. 0) THEN
               WRITE (*,*) ' Density difference requested',
     +                     ' but no root number specified!'
               call errex
            ELSE IF (iLst .NE. 2) THEN
               WRITE (*,*) ' Incorrect root number specified',
     +                     ' for density difference!'
               WRITE (*,*) ' Must be in form XX (i.e. 01, 02, etc)!'
               call errex
            ELSE
               DLABEL="TDENSITY"
               iroot = args(3)(1:iLst)
            ENDIF
          ENDIF
        ELSE
          DLABEL = 'TDENSITY'
        ENDIF

CSSS        WRITE (*,*) 'DLABEL - ', DLABEL, LNAT_ORBS
        IF ( .NOT. LNAT_ORBS) THEN
           iener=inext
           iocc=iener+nmo*iintfp
           iorb=iocc+nmo+mod(nmo,2)
           iorbr=iorb+nao*nmo*iintfp
           imom =iorbr+nao*iintfp
           inext=imom+nao
           
           if(inext-i0.gt.maxcor)call insmem('rdorb-f',inext-i0,maxcor)

           call molden_rdorb(icore(iener),icore(iocc),icore(iorb),
     &          icore(iorbr),icore(imom),icore(inext),nao,nmo,
     &          maxcor-(inext-i0),iuhf,iexx,iunit)

        ELSE
           iener=inext
           iocc=iener+nmo*iintfp
           iorb=iocc+nmo+mod(nmo,2)
           igden=iorb+nao*nmo*iintfp
           ivec=igden+nmo*nmo*iintfp
           ieden=ivec+nmo*nmo*iintfp
           inext=ieden+nmo*nmo*iintfp

           if(inext-i0.gt.maxcor)call insmem('nlorb ',inext-i0,maxcor)

           call molden_nlorb (icore(iener),icore(iocc),icore(iorb),
     &       icore(igden),icore(ivec),icore(ieden),
     &       icore(inext),nao,nmo,maxcor-(inext-i0),iuhf,iexx,iunit,
     &       DLABEL,iroot,dens_diff)
        ENDIF
      endif
c-------------------------------------------------------------------------
c vibrational frequencies
c
      if(iflags(54).ne.0) then
        write(iunit,*)  
        write(iunit,'(a)')'[FREQ]'

        call getrec(-1,'JOBARC','NREALATM',1,nreal_atoms)
        nvib=3*nreal_atoms
 
        ifreq=inext
        ifreqco=ifreq+nvib*iintfp
        inormmd=ifreqco+3*nreal_atoms*iintfp
        inext=inormmd+nvib*3*nreal_atoms*iintfp

        if(inext-i0.gt.maxcor)call insmem('rdvib-f',inext-i0,maxcor)
        call A2rd_vib(nreal_atoms,icore(iatchrg),icore(ifreq),
     &    icore(icoord),icore(inormmd),Nimag,nvib,iunit,.true.)

      endif

c ----------------------------------------------------------------------

      if (iunit.ne.6) then
         close(iunit)
         print *, '@MOLDEN_MAIN: successfully created ',szFile
      end if

      return
      end

