      subroutine A3_ecp_basis_init_4grd(Natoms, Namat, Iqmstr,
     &                                  Jcostr, Nucstr, Nrcstr, Jstrt, 
     &                                  Jrs, Expa, Nalpha, Cont, Npcoef,
     &                                  Nharm)
c
c-----------------------------------------------------------------------      
c     subroutine ECP_basis_init fills the common-blocks for subroutine 
c     version for general contracted basissets
C
c     input:
c     namstr: names of the atoms (strange storage as integer!)
c     iqmstr: max. l-quantum-number+1 of a shell
c     jcostr: number of shells
c     nucstr: number of the uncontracted functions of the shells > NUCO
c             (exponents)
c     nrcstr: number of blocks(general contracted) > NRCO
c             (contractions)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      integer eor,and,or,dstrt
C-----------------------------------------------------------------------
#include "ecp.par"
#include "ecp.com"
C-----------------------------------------------------------------------
      character*4 namat(Max_centers)
      logical nharm

      Dimension iqmstr(Max_centers),jcostr(Max_centers,maxjco),
     &          nucstr(Max_centers,Maxang,Maxjco),
     &          nrcstr(Max_centers,Maxang,Maxjco),
     &          Jstrt(Max_shells),Jrs(Max_shells),Expa(Nalpha),
     &          Cont(Npcoef)
C
      Dimension
     &          ipeoff(Max_centers,Maxang,Maxjco),
     &          iccoff(Max_centers,Maxang,Maxjco),
     &          newoff(Max_centers,Maxang,Maxjco),
     &          numpri(Max_centers,Maxang,Maxjco),
     &          imin(Maxang),imax(Maxang),nfcts(Maxang),
     &          ifeld(Max_centers)
C
c-----------------------------------------------------------------------
c     gout                : intbatch
c     ipeoff              : offset for primitive exponents
c     iccoff              : offset for contraction coefficients
c     newoff              : offsets for primitives
c     numpri              : number of primitive functions of a shell
c     imin,imax           : lower(upper) offset of a cartesian component
c     nfcts               : number of types of a cartesian component 
c     ifeld               : checkarray 
c-----------------------------------------------------------------------
CSSS      common /infoa / natoms
CSSS      common /symshe/ mulsh(Mxtnsh)
CSSS      common /modez / zetm(Mxtnsh)
CSSS      common /pairij/ ipq(ndi4+1)
CSSS      common /typinf/ ltmax,ltdrv,nftnft,nftdrv
CSSS      common /powers/ jx(nftmax),jy(nftmax),jz(nftmax),
CSSS     1                ix(nftmax),iy(nftmax),iz(nftmax)
C------------------------------------------------------------------------
c Create offsets of primitive exponents and contraction coefficients
C The NMax refers 
C
      icount=0
      call izero(ifeld,Mxatms)
CSSS      call izero(ipeoff,Mxatms*Maxang*Mxshel)  
CSSS      call izero(iccoff,Mxatms*Maxang*Mxshel)  

      do  iat=1,Natoms
        do  ilq=1,iqmstr(iat)
          do  ijco=1,jcostr(iat,ilq)

              icount=icount+1
C
              if (ifeld(iat).ne.1) then

                 ipeoff(iat,ilq,ijco)=jstrt(icount)
                 iccoff(iat,ilq,ijco)=nstrt(icount)
#ifdef _DEBUG_LVL0
      write(*,1) iat,ilq,ijco,ipeoff(iat,ilq,ijco)
   1  format('iat: ',i3,'ilq: ',i3,' ijco: ',i3,' ipeoff: ',i3)
      write(*,2) iat,ilq,ijco,iccoff(iat,ilq,ijco)
   2  format('iat: ',i3,'ilq: ',i3,' ijco: ',i3,' iccoff: ',i3)
#endif
              endif
          enddo
        enddo 
C
        ifeld(iat)=1 
C
      enddo

C create shell-structure and arrays expnt and contr

CSSS      call zero(expnt,Mxtnpr)
CSSS      call zero(contr,Mxtnpr*Mxtnpr)
CSSS      call izero(numpri,Mxatms*Maxang*Maxprim)

      icount=0
      kcount=0
      ncount=0

      do iat=1,Natoms
        do ilq=1,iqmstr(iat)
          do  ijco=1,jcostr(iat,ilq)
C
            ncount=ncount+1

            do inuco=1,nucstr(iat,ilq,ijco)

              icount=icount+1

              do inrco=1,nrcstr(iat,ilq,ijco)
C
                 contr(icount,inrco)=CONT(ipeoff(iat,ilq,ijco)+inuco+
     &                               nucstr(iat,ilq,ijco)*(inrco-1))

#ifdef _DEBUG_LVL0
      write(*,3) icount,inrco,contr(icount,inrco)
   3  format('icount: ',i3,' inrco: ',i3,' contr: ',f12.8)
#endif
              enddo

              expnt(icount)=EXPA(ipeoff(iat,ilq,ijco)+inuco)

#ifdef _DEBUG_LVL0
      write(*,4) icount,expnt(icount)
   4  format('icount: ',i3,' expnt: ',f12.8)
#endif
              kcount=kcount+1

            enddo

c           numpri(iat,ilq,ijco)=numpri(iat,ilq,ijco)+kcount
            numpri(iat,ilq,ijco)=kcount

#ifdef _DEBUG_LVL0
      write(*,5) iat,ilq,ijco,numpri(iat,ilq,ijco)
   5  format('iat: ',i3,' ilq: ',i3,' ijco: ',i3,' numpri: ',i3)
#endif
            kcount=0
            numcon(ncount)=nrcstr(iat,ilq,ijco)

#ifdef _DEBUG_LVL0
      write(*,6) ncount,numcon(ncount)
   6  format('ncount: ',i3,' numcon: ',i3)
#endif
          enddo

c         numpri(iat,ilq)=numpri(iat,ilq)+kcount
c         write(*,5) iat,ilq,numpri(iat,ilq)
c5        format('iat: ',i3,' ilq: ',i3,' numpri: ',i3)
c         lcount=lcount+1
c         kcount=0

        enddo

c       lcount=0

      enddo
C
c create new offsets of primitive exponents
C
CSSS      call izero(newoff,Mxatms*Maxang,Maxprim)

      khelp=1

      do iat=1,NMAX
        do  ilq=1,iqmstr(iat)
          do ijco=1,jcostr(iat,ilq)

            newoff(iat,ilq,ijco)=khelp

#ifdef _DEBUG_LVL0
      write(*,99) iat,ilq,ijco,newoff(iat,ilq,ijco)
  99  format('iat: ',i3,' ilq: ',i3,' ijco: ',i3,' newoff: ',i3)
#endif
            khelp=khelp+numpri(iat,ilq,ijco)

          enddo
        enddo
      enddo

c-----------------------------------------------------------------------
c     fill in common infoa
c-----------------------------------------------------------------------

CSSS      call zero(xyz,3*Mxatms)

      natoms=NMAX

CSSS      do  i=1,natoms
CSSS        do  j=1,3
CSSS          xyz(j,i)=CORD(i,j)
CSSS        enddo
CSSS      enddo

CSSS      do  i=1,natoms
CSSS        charg(i)=CHARGE(i)
CSSS      enddo
c-----------------------------------------------------------------------
c     fill in common /typinf/
c-----------------------------------------------------------------------

      ltmax=1
      do 10 iat=1,NMAX
        if (ltmax.le.iqmstr(iat)) then 
            ltmax=iqmstr(iat)
        endif  
   10 continue      

      nftdrv=ltmax*(ltmax+1)*(ltmax+2)/6

c-----------------------------------------------------------------------
c     fill in common /nshell/
c-----------------------------------------------------------------------

c Doing the k's from common nshel
      
      if (NHARM) then

c spherical harmonics are used!

        do l=1,Maxang
           imin(l)=(l-1)*l*(l+1)/6+1
           imax(l)=l*(l+1)*(l+2)/6 
           nfcts(l)=2*l-1
        enddo

      else 

        do l=1,Maxang
           imin(l)=(l-1)*l*(l+1)/6+1
           imax(l)=l*(l+1)*(l+2)/6 
           nfcts(l)=l*(l+1)/2
        enddo
      endif

      nbf=0
      nshll=0
      indsh=0
C
c loop over atoms

      do iat=1,NMAX

c loop over l-quantum-numbers

        do ilq=1,iqmstr(iat)

c loop over contracted basisfunctions

          do ijco=1,jcostr(iat,ilq)

             indsh=indsh+1
             mulilq=ISTBAO(indsh)
             nshll=nshll+1
             katom(nshll)=iat
             ktype(nshll)=ilq
             kmini(nshll)=imin(ilq)
             kmaxi(nshll)=imax(ilq)
             ndegen(nshll)=nfcts(ilq)             
             kbfn(nshll)=nbf+1
             nprims(nshll)=numpri(iat,ilq,ijco)
             kprim(nshll)=newoff(iat,ilq,ijco)
             nbf=nbf+ndegen(nshll)*fmult(mulilq)

          Enddo
        Enddo
      Enddo

      nshell=nshll

#ifdef _DEBUG_LVL0
      write(*,*)
      do  i=1,nshell
         write(*,*) 'nshell=',i
         write(*,3913) i,katom(i),i,ktype(i)
3913     format('katom<',i2,'>  ',i2,' ktype<',i2,'>  ',i2)
         write(*,3914) i,kprim(i),i,kbfn(i)
3914     format('kprim<',i2,'>  ',i2,' kbfn<',i2,'>   ',i2)
         write(*,3915) i,nprims(i),i,ndegen(i)
3915     format('nprims<',i2,'> ',i2,' ndegen<',i2,'> ',i2)
         write(*,3916) i,kmini(i),i,kmaxi(i)
3916     format('kmini<',i2,'>  ',i2,' kmaxi<',i2,'>  ',i2)
         write(*,3917) i, mulsh(i)
3917     format('mulsh<',i2,'>  ',i2)
         write(*,*)
      enddo
      write(*,*) 'nbf=',nbf
      write(*,*)
#endif

c-----------------------------------------------------------------------
c     fill in common modez
c-----------------------------------------------------------------------

CSSS      do ii=1,nshell
CSSS        zetm(ii)=expnt(kprim(ii))
CSSS        do ipr=kprim(ii)+1,kprim(ii)+nprims(ii)-1
CSSS          zetm(ii)=min(expnt(ipr),zetm(ii))
CSSS        enddo
CSSS      enddo
c-----------------------------------------------------------------------
c     fill in common pairij
c-----------------------------------------------------------------------

CSSS      call setipq(ipq,ndi4)

c-----------------------------------------------------------------------
c     fill in common pairij
c-----------------------------------------------------------------------

CSSS      call monom(ltmax,nft,jx,jy,jz,iadr,nftmax,1,0)

c-----------------------------------------------------------------------
c     fill common /palcal/
c-----------------------------------------------------------------------
CSS      nshdim=nshell
c      TLA=1.0d-12
c      toll=-log10(TLA)+8.0
CSSS       toll=25.0    
C
      return
      end 
