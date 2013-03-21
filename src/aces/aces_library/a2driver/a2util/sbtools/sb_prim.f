
c This routine reads the MOL file and initializes data in /MOL/.
c It should be called only by sb_com_mol(), but its main purpose is
c to recast icore() and dcore() addresses into convenient variables.
c
c Loads primitive coefficients of the basis.

      subroutine sb_prim(comppopv,compmemb,nshellatom,nshellprim,
     &                   nshellao,primoff,aooff,alpha,primcoef,
     &          dpmcoef,nprimao,nprimaol)
      implicit none



c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)




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



      integer comppopv(*),compmemb(*),nshellatom(*),
     &        nshellprim(maxshell,*),nshellao(maxshell,natoms),
     &        primoff(maxshell,*),aooff(maxshell,*)
      double precision alpha(*),primcoef(totprim,*)
       integer iip,ikk,ill,iloo,ioo,ijj,
     &  nprimao(naobasfn),nprimaol(naobasfn),incre

      double precision nthres,
     &   dpmcoef(totprim,naobasfn)
      double precision pi,s,tmp,xnorm,ai,aj,dpow
      integer
     &    iorbit,ishell,ishelloff,ishellang,isub,iang,
     &    iequiv,nsubshell(maxangshell),nshell,iatom,iprimoff,iaooff,
     &    isubprim,isubao,i,nline,iprim0,iprim1,iprim,iao0,iao1,iao,
     &    jprim,jprim0,jao0,jprim1,jao1,jao,icount,iequivatom,iunit
      character*80 line,molfil

      call callstack_push('SB_PRIM')

      pi = acos(-1.d0)
      nthres=10.d0**(-12.d0)
 9000 format(a25,12i5)
 9010 format(2i5)
 9020 format(4f18.10)

      iunit=3
      call gfname('MOL',molfil,i)
      open(unit=iunit,file=molfil,form='formatted',status='old')

c Set nshellprim and nshellao
      rewind(iunit)
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

c iorbit     the orbit number
c icount     the pointer to compmemb for the current atom of the current orbit
c iatom      the atom at that location in compmemb
c iequivatom the current equivalent atom

      icount=1
      do iorbit=1,compnorb

         iatom=compmemb(icount)
         icount=icount+1
         read(iunit,9000) line,nshell,
     &                    (nsubshell(ishell),ishell=1,nshell)

         do iip=1,nshell
         end do
         read(iunit,'(a)') line
         ishellang=0
         ishelloff=1
         do ishell=1,nshell
c ishelloff=1,2,5,11,21,...  (where the first shell of ishellang starts)
c ishellang=1,3,6,10,15,...  (1s, 3p, 6d, 10f, etc)
            ishelloff=ishelloff+ishellang
       
            ishellang=ishell*(ishell+1)/2
             
     
            do isub=1,nsubshell(ishell)
               read(iunit,9010) isubprim,isubao
                 

               do iang=ishelloff,ishelloff+ishellang-1
                  nshellprim(iang,iatom)=nshellprim(iang,iatom)+isubprim

                  nshellao(iang,iatom)=nshellao(iang,iatom)+isubao

               end do
               nline=(isubao-3)/4
                
               if ((isubao-3).gt.(nline*4)) nline=nline+1
               nline=(nline+1)*isubprim
               do i=1,nline
                  read(iunit,'(a)') line
               end do
            end do
         end do

c Copy primitive and ao info to equivalent atoms
         if (comppopv(iorbit).gt.1) then
            do iequiv=2,comppopv(iorbit)
               iequivatom=compmemb(icount)
               icount=icount+1
               do ishell=1,nshellatom(iatom)
                  nshellprim(ishell,iequivatom)=nshellprim(ishell,iatom)
                  nshellao(ishell,iequivatom)=nshellao(ishell,iatom)
               end do
            end do
         end if

c     end do iorbit=1,compnorb
      end do

c Set up primoff and aooff
      iprimoff=1
      iaooff=1
      do iatom=1,natoms
         do ishell=1,nshellatom(iatom)
            primoff(ishell,iatom)=iprimoff
            iprimoff=iprimoff+nshellprim(ishell,iatom)
            aooff(ishell,iatom)=iaooff
            iaooff=iaooff+nshellao(ishell,iatom)
         end do
      end do

c Read in the exponents and primitive coefficients (renormalizing them)
      rewind(iunit)
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

      icount=1
      do iorbit=1,compnorb
         iatom=compmemb(icount)
         read(iunit,9000) line,nshell,
     &                    (nsubshell(ishell),ishell=1,nshell)
         read(iunit,'(a)') line
         ishellang=0
         ishelloff=1
         do ishell=1,nshell
            ishelloff=ishelloff+ishellang
            ishellang=ishell*(ishell+1)/2
c iprim0,iprim1 : the first and last primitives in the current set
c iao0,iao1     : the first and last AOs in the current set
            iprim0=primoff(ishelloff,iatom)
            iao0=aooff(ishelloff,iatom)
            do isub=1,nsubshell(ishell)
               read(iunit,9010) isubprim,isubao
               iprim1=iprim0+isubprim-1
               iao1=iao0+isubao-1
c Read all exponents and coefficients

               do iprim=iprim0,iprim1
                  read(iunit,9020) alpha(iprim),
     &                             (primcoef(iprim,iao),iao=iao0,iao1)
                 do iao=iao0,iao1
               end do
               end do
c Renormalize

               do iao=iao0,iao1
                  s=0.d0
                  dpow=dble(ishell)+0.5d0
                  do iprim=iprim0,iprim1
                     ai=alpha(iprim)
                     do jprim=1,iprim
                        aj=alpha(jprim)
                        tmp=primcoef(iprim,iao)
     &                     *primcoef(jprim,iao)
     &                     *(2.d0*sqrt(ai*aj)/(ai+aj))**dpow
                        s=s+tmp
                        if (iprim.ne.jprim) s=s+tmp
                     end do
                  end do
                  xnorm=((0.5d0/pi)**(3.d0*0.25d0))/sqrt(s)
                  dpow=0.5d0*dpow
                  do iprim=iprim0,iprim1
                     primcoef(iprim,iao)=primcoef(iprim,iao)
     &                                  *xnorm
     &                                  *(4.d0*alpha(iprim))**dpow
                  end do
               end do
               iao0=iao0+isubao
               iprim0=iprim0+isubprim
            end do
         end do

c Copy to equivalent blocks (given iorbit,ishell)
c This time we have to loop over ALL equivalent atoms since we'll also
c need to fill in the remaining shells of the first atom.
         do iequiv=1,comppopv(iorbit)
            iequivatom=compmemb(icount)
            icount=icount+1
            ishellang=0
            ishelloff=1
            do ishell=1,nshell
               ishelloff=ishelloff+ishellang
               ishellang=ishell*(ishell+1)/2
               iprim0=primoff(ishelloff,iatom)
               iao0=aooff(ishelloff,iatom)
               do iang=ishelloff,ishelloff+ishellang-1
                  jprim0=primoff(iang,iequivatom)
                  jprim1=jprim0+nshellprim(iang,iequivatom)-1
                  jao0=aooff(iang,iequivatom)
                  jao1=jao0+nshellao(iang,iequivatom)-1
                  do jprim=jprim0,jprim1
                     iprim=(jprim-jprim0)+iprim0
                     alpha(jprim)=alpha(iprim)
                     do jao=jao0,jao1
                        iao=(jao-jao0)+iao0
                        primcoef(jprim,jao)=primcoef(iprim,iao)
                     end do
c                 end do jprim=jprim0,jprim1
                  end do
c              end do iang=ishelloff,ishelloff+ishellang-1
               end do
c           end do ishell=1,nshell
            end do
c        end do iequiv=1,comppopv(iorbit)
         end do
c     end do iorbit=1,compnorb
      end do

      close(iunit)

c          iloo=1
c          nprimao(1)=1
c           nprimaol(1)=0 
        do ikk=1,naobasfn
           incre=0
        do ill=1,totprim
           if(abs(primcoef(ill,ikk)) .gt. nthres) then
               if(incre .eq.0) then
                   nprimao(ikk)=ill
                   nprimaol(ikk)=ill
                   incre=incre+1
                else
                   nprimaol(ikk)=ill
                end if
            end if
         end do
         end do
                  













c              if(ill .le. nprimaol(ikk-1)) then
                
c                 if(incre .eq.0) then
c                    nprimao(ikk)=ill
c                  end if
c                  incre=incre+1 
c                  iloo=ill
c               dpmcoef(iloo,ikk)=primcoef(ill,ikk)
c              iloo=iloo+1
c              else
c                dpmcoef(iloo,ikk)=primcoef(ill,ikk)
c                iloo=iloo+1
c              end if
c          end if
c        end do
c             nprimaol(ikk)=iloo-1
c             write(*,*) ikk,nprimaol(ikk) 
c             if(ikk .lt. naobasfn) then
c               nprimao(ikk+1)=iloo
c            end if
c        end do

c       call putrec(1,'JOBRAC','CCOEFFC',totprim*naobasfn*iintfp
c     &             ,dpmcoef)
     
       call putrec(1,'JOBARC','NPRIMST',naobasfn,nprimao)
       call putrec(1,'JOBARC','NPRIMLT',naobasfn,nprimaol) 
      
      call callstack_pop
      return
      end

