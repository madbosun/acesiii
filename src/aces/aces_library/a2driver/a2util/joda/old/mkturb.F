      subroutine mkturb(basnam,geom,zsym,atmass,iatnum,natms,
     &           nontyp,ispher)
c
c --- <mkturb> prepares the shellinformation (i.g. common /nshel / and /infao /)
c              of the TURBOMOLE interface 
cend
      implicit double precision (a-h,o-z)

#include "mxatms.par"
      parameter (ndi13  = 1000)
      parameter (ndi9   = 250)
      parameter (ndi10  = mxatms)
      parameter (ndilmx = 5)
c
      character*8 atsymb(ndi10)
      character*5 tmfil
      character*(*) basnam(mxatms)
csb 1/97 Hold long basis set names
      integer baslen
      parameter (baslen=80)
      CHARACTER*(baslen) TITEL(ndi10)

      character*5 zsym(natms)
      character*4 fpgrp, bpgrp, pgrp
c
      integer shellang,shellloc,shellsiz,shellprm 
      integer turbo
      parameter (turbo  = 66)
      parameter (tmfil  = 'TURBO')
c
      logical opn
c
      dimension geom(3*mxatms),atmass(natms)
      dimension iatnum(natms)
      dimension imin(ndilmx),imax(ndilmx),nfcts(ndilmx),nfctst(ndilmx)
      dimension ktypef(ndi9),katomf(ndi9),kprimf(ndi9),nprimf(ndi9)
c      dimension expntf(ndi13),contrf(ndi13)
      dimension nmtshl(ndi10)
      dimension xyz(3,ndi10),charge(ndi10),wmass(ndi10)
      dimension expnt(ndi13),contr(ndi13),katom(ndi9),ktype(ndi9),
     &          kprim(ndi9),kbfn(ndi9),kmin(ndi9),kmax(ndi9),
     &          nprims(ndi9),ndegen(ndi9),zetm(ndi9)
      dimension krbfn(ndi9)
      dimension npople(5)
      dimension ibasex(mxatms)
c
      common /syminf/ fpgrp,bpgrp,pgrp
      common /tmpinf/ nshlat(ndi10),iatloc(ndi10),titel,
     &                iofshl(ndi10),iofprm(ndi9)
      common /molstr1/ shellang(100),shellloc(100),shellsiz(100),
     &                 shellprm(100),basisexp(10000),basiscnt(10000)
c
c --- set imin,imax arrays which contain start and end indices of
c     canonical ordered cartesian gaussian basis functions up to mxltyp
c     { example : imin(s)=imax(s)=1 ; imin(d)=5,imax(d)=10 ... }
c     set nfcts array, which contain number of (reducible) cartesian
c     gaussian functions of a given l quantum number
c
      if (ispher.eq.1) then
        npople(1) = 1
        npople(2) = 3
        npople(3) = 5
        npople(4) = 7
        npople(5) = 9
      elseif (ispher.eq.0) then
        npople(1) = 1
        npople(2) = 3
        npople(3) = 6
        npople(4) = 10
        npople(5) = 15
      endif
      do 100 l=1,ndilmx
        imin(l)=(l-1)*l*(l+1)/6+1
        imax(l)=l*(l+1)*(l+2)/6
        nfcts(l)=l*(l+1)/2
        nfctst(l)=2*l-1
  100 continue
c
c --- init ibasex
c
      call izero(ibasex,natms)
c
c --- get the information needed by shloop and paloop 
c
      icntat=0
      nshll=0
      iprtot=0
      mprim=1
c
c --- loop over ALL atoms
c
      ibold=0
      do 1000 ia=1,natms
        if (index(basnam(ia),':').ne.0) then
          icntat=icntat+1
          atsymb(icntat)=zsym(ia)
          wmass(icntat)=atmass(ia)
          charge(icntat)=float(iatnum(ia))
          icount=0
          iofst=3*(ia-1)
          do 1100 jj=1,3
            icount=icount+1
            xyz(jj,icntat)=geom(iofst+icount)
 1100     continue
          do 1200 iats=1,nontyp
            if (basnam(ia).eq.titel(iats)) then
              nmshl=nshlat(iats)
              nmatm=iatloc(iats)
              ibasex(icntat)=iats
            endif
 1200     continue
          icke=0
c
c ---     loop over the shells of the atom
c
          if (ibasex(icntat).gt.ibold) then
            ibold=ibasex(icntat)
            do 2000 is=1,nmshl
              ioffs=iofshl(nmatm)
              iproff=iofprm(ioffs+is)
              icke=icke+shellsiz(ioffs+is)
c              nmtshl(icntat)=icke
              nmtshl(ibasex(icntat))=icke
              ncont=shellsiz(ioffs+is)
              nprim=shellprm(ioffs+is)
              iprshl=0
c
c ---       loop over the blocks of contraction coefficients
c
              do 3000 ic=1,ncont
c
c ---         loop over the primitives
c
                do 4000 ip=1,nprim
                  if (abs(basiscnt(iproff+(ic-1)*nprim+ip)).
     &            gt.1.0d-15) then
                    iprshl=iprshl+1
                    iprtot=iprtot+1
                    contr(iprtot)=basiscnt(iproff+(ic-1)*nprim+ip)
                    expnt(iprtot)=basisexp(iproff+(ic-1)*nprim+ip)
                  endif
 4000           continue
                nshll=nshll+1
                nprims(nshll)=iprshl
                kprimf(nshll)=mprim
                mprim=mprim+nprims(nshll)
                katom(nshll)=icntat
                ktype(nshll)=shellang(ioffs+is)+1
                iprshl=0
 3000         continue
 2000       continue
          endif
        endif
 1000 continue
c      natom=icntat
c
c --- find total number of primitives and atoms
c
      natom=icntat
      numprn=iprtot    
c
c --- find the highest quantum number of the shells
c
      ltmax=1
      do 200 i=1,nshll
        ltmax=max(ktype(i),ltmax)
  200 continue

      icnt=0
      nbf=0 
      nsymbf=0 
      ipoff=0
      mprim=1
c
c --- sort the shells by the angular momentum numbers
c
c --- loop over the angular momentum number
c
      do 5000 il=1,ltmax
c
c ---   loop over atoms
c
        do 5500 iato=1,natom
          icons0=0
          do 5550 ibas=1,ibasex(iato)-1
 5550      icons0=icons0+nmtshl(ibas)
c
c ---     loop over the shells
c
c          do 6000 is=1,nshll
          do 6000 is=1,nmtshl(ibasex(iato))
            icons=icons0+is
            if (ktype(icons).eq.il) then
              icnt=icnt+1
              ktypef(icnt)=ktype(icons)
              katomf(icnt)=iato
              ndegen(icnt)=nfcts(ktype(icons))
              kmax(icnt)=imax(ktype(icons))
              kmin(icnt)=imin(ktype(icons))
              kbfn(icnt)=nbf+1
              krbfn(icnt)=nsymbf
              nbf=nbf+ndegen(icnt)
c            nsymbf=nsymbf+nfctst(ktype(icons))
              nsymbf=nsymbf+npople(ktype(icons))
              kprim(icnt)=kprimf(icons)
              nprimf(icnt)=nprims(icons)
              mprim=mprim+nprimf(icnt)
c
c ---       loop over the primitives of the shell
c
c            do 6500 ip=1,nprims(is)
c              expntf(ip+ipoff)=expnt(kprim(is)+ip-1)
c              contrf(ip+ipoff)=contr(kprim(is)+ip-1)
c 6500       continue
c            ipoff=ipoff+nprims(is)
            endif
 6000     continue
 5500   continue
 5000 continue
c
c --- determine minimum exponents of contracted basis function blocks
c
      nshll=icnt
      do 300 ii=1,nshll
        zetm(ii)=expnt(kprim(ii))
        do 400 ipr=kprim(ii)+1,kprim(ii)+nprimf(ii)-1
          zetm(ii)=min(expnt(ipr),zetm(ii))
  400   continue
  300 continue
c
c --- store nbf to the JOBARC-file for the VTRAN-module
c
      call putrec(10,'JOBARC','NBFINCAO',1,nbf)
c
c --- now we are ready to store the information on file
c
c --- open the ACES2-TURBOMOLE interface file
c
      inquire(file=tmfil,exist=opn)
      if (opn) then
        open(file=tmfil,unit=turbo,form='formatted',status='old')
        close(unit=turbo,status='delete')
      endif
      open(file=tmfil,unit=turbo,form='formatted',status='new')
c
      write(turbo,'(a)') 'TURBOACES'
c
c --- 1. information about the atoms
c
      write(turbo,'(a)') 'ATOMINFO'
      write(turbo,'(i5)') natom
      do 7000 ia=1,natom
        write(turbo,7010) atsymb(ia),charge(ia),wmass(ia),
     &                    nmtshl(ibasex(ia))
 7010   format(a8,f7.3,f20.15,i4)
        write(turbo,7020) (xyz(i,ia),i=1,3)
 7020   format(3f20.13)
 7000 continue
c
c --- 2. information about the symmetry of the molecule
c
      write(turbo,'(a)') 'SYMMETRYINFO'
      write(turbo,'(3a4)') fpgrp,bpgrp,pgrp
c
c --- 3. information about the shells
c
      write(turbo,'(a)') 'SHELLINFO'
      write(turbo,'(4i5)') nbf,nsymbf,nshll,ltmax
      do 8000 is=1,nshll
        write(turbo,8010) katomf(is),ktypef(is),kprim(is),nprimf(is),
     &                    kbfn(is),krbfn(is),ndegen(is),kmin(is),
     &                    kmax(is),zetm(is)
 8010   format(9i5,f25.15)
 8000 continue
c
c --- 4. information about the primitives
c
      write(turbo,'(a)') 'PRIMINFO'
c      write(turbo,'(i5)') ipoff
      write(turbo,'(i5)') numprn
      do 9000 ip=1,numprn
        write(turbo,9010) expnt(ip),contr(ip)
 9010   format(2f25.15)
 9000 continue
c
      write(turbo,'(a)') 'FINISH'
      close(unit=turbo,status='keep')
c
      return
      end
