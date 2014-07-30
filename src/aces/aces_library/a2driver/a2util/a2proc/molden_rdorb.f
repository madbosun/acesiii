      subroutine molden_rdorb(ener,iocc,orb,orb_r,iang,scr,nao,nmo,
     &                        maxcor,iuhf,iexx,iunit)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision ener(nmo),orb(nao,nmo),scr(maxcor),
     &                 orb_r(nao)
      character*8 cscfener(2)
      character*8 cscforb(2)
      character*8 cexxcoef(2)
      character*8 cnumdrop(2)
      character*5 sptype(2)
      integer iocc(nmo), iang(nao)
      integer  idrppop(8),idrpvrt(8)



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



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
      parameter (one=1.0D0)
      parameter (zilch=0.0D0)
      data cscfener /'SCFEVLA0','SCFEVLB0'/
      data cscforb  /'SCFEVCA0','SCFEVCB0'/
      data cexxcoef /'EXXCOEFA','EXXCOEFB'/
      data cnumdrop /'NUMDROPA','NUMDROPB'/  
      data sptype   /'Alpha','Beta '/
      call aces_com_info
      call aces_com_syminf
      call aces_com_sym
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c the mo ordering is important for getting the occupations right
c
c depending on the place in the program (if the user is not using xaces2)
c it may be either scf or correlated
c
c    so            mo                        mo
c ao CMP2ZMAT   so cscforb(ispin)    =    ao AOBASMOS
c
c the so ordering is zmat
c the mo ordering is correlated if calc is greater than scf 
c or if this is a vibrational calculation
c
      i000=1
      i010=i000+nao*nmo
      i020=i010+nmo**2

      if(i020.gt.maxcor/iintfp)
     &  call insmem('RDORB-F',i020*iintfp,maxcor)

      call getrec(-1,'JOBARC','CMP2ZMAT',nao*nmo*iintfp,scr(i000))
c      if(iexx.eq.1) then
c          write(*,*)" number of AOs: ",nao
c          write(*,*)" number of MOs: ",nmo
c         write(*,*)" cmp2zmat matrix is: "
c         do i=1,nao*nmo
c           write(*,*)scr(i000+i-1)
c         enddo
c      endif
      do 10 ispin=1,iuhf+1
        if(iexx.eq.1)then
           call getrec(-1,'JOBARC',cexxcoef(ispin),nao*iintfp,
     $          scr(i010))
c           write(*,*)" the exx potential coefficents are: "
c           do i=1,nao
c             write(*,*) scr(i010+i-1)
c           enddo
c           do i=2,nmo
c              call dcopy(nao,scr(i010),1,scr(i010+(i-1)*nao),1)
c           enddo
        else
           call getrec(20,'JOBARC',cscforb(ispin),
     $          nmo**2*iintfp,scr(i010))

        endif

        call xgemm('n','n',nao,nmo,nmo,one,scr(i000),nao,scr(i010),
     &    nmo,zilch,orb,nao)

        call getrec(-1,'JOBARC',cscfener,nmo*iintfp,ener)
c
c fill occupation numbers
c
cYAU - iflags(52) corresponds to geometry search step size controls...
        if((iflags(2).gt.0).or.(iflags(52).gt.0))then
c
c correlated ordering
c if dropmo change to full dimension
c
          call getrec(20,'JOBARC',cnumdrop(ispin),1,idrop)
          if(idrop.gt.0) then
            call getrec(-1,'JOBARC','NDROPPOP',nirrep,idrppop)
            call getrec(-1,'JOBARC','NDROPVRT',nirrep,idrpvrt)
            do 20 irrep=1,nirrep
              pop(irrep,ispin)=pop(irrep,ispin)+idrppop(irrep)
              nocco(ispin)=nocco(ispin)+idrppop(irrep)              
              vrt(irrep,ispin)=vrt(irrep,ispin)+idrpvrt(irrep)
              nvrto(ispin)=nvrto(ispin)+idrpvrt(irrep)
20          continue
          endif
          if(iuhf.eq.0)then
            do 30 imo=1,nocco(ispin)
              iocc(imo)=2
30          continue
          else
            do 35 imo=1,nocco(ispin)
              iocc(imo)=1
35          continue
          endif
          do 40 imo=nocco(ispin)+1,nmo
            iocc(imo)=0
40        continue 
        else
c
c scf ordering dropmo 
c since vtran hasn't run dropmo should be irrelevant here
c
          ioff=1
          if(iuhf.eq.0)then
            do 50 irrep=1,nirrep
              do 51 imo=ioff,ioff+pop(irrep,ispin)-1
                iocc(imo)=2
51            continue
              ioff=ioff+pop(irrep,ispin)
              do 52 imo=ioff,ioff+vrt(irrep,ispin)-1
                iocc(imo)=0
52            continue
              ioff=ioff+vrt(irrep,ispin)
50          continue
          else
            do 60 irrep=1,nirrep
              do 61 imo=ioff,ioff+pop(irrep,ispin)-1
                iocc(imo)=1
61            continue
              ioff=ioff+pop(irrep,ispin)
              do 62 imo=ioff,ioff+vrt(irrep,ispin)-1
                iocc(imo)=0
62            continue
              ioff=ioff+vrt(irrep,ispin)
60          continue
          endif
        endif
C
        do 80 imo=1,nmo
          write(iunit,90)' Ene=   ',ener(imo)
90        format(A,F8.4)
          write(iunit,100)' Spin= ' // sptype(ispin)
100       format(A)
          write(iunit,110)'  Occup=   ',iocc(imo),'.000000'
110       format(A,I1,A)
C
C Reorder the F functions (perhaps g too ..), Peter Szalay identified
C the problem and proivided the reorder routines. 
C Ajith Perera, 06/2011.
C
          call getrec(-1, 'JOBARC', 'ANMOMBF0', nao, iang)
          call reorder(orb(1, imo), orb_r, iang, nao)
          do 120 iao=1,nao
CSSS             write(iunit,130)iao,'  ',orb(iao,imo)
            write(iunit,130)iao,'  ',orb_r(iao)
130         format(I4,A,F12.6)
120       continue
80      continue 
10    continue
      return
      end




