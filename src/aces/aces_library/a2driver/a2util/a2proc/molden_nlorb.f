      subroutine molden_nlorb(ener,iocc,orb,dens,nlorb,edens,scr,
     &                        nao,nbas,maxcor,iuhf,iexx,iunit,
     +                        DLABEL,iroot,dens_diff)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision ener(nbas),orb(nao,nbas),scr(maxcor)
      double precision dens(nbas,nbas),edens(nbas,nbas),nlorb(nbas,nbas)
      double precision nelec
      logical  dens_diff
      character*8 DLABEL
      character*2 iroot
      character*8 cscfener(2)
      character*8 cscforb(2)
      character*8 cexxcoef(2)
      character*8 cnumdrop(2)
      character*5 sptype(2)
      double precision  iocc(nbas)
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
      parameter (DENS_THRESH=1.0D-08)
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
      i010=i000+nbas*nbas
      i020=i010+nbas*nbas
      i030=i020+nbas*nbas
      iend=i030+2*nbas*nbas

      if(i030.gt.maxcor/iintfp)
     &  call insmem('NLORB-F',iend*iintfp,maxcor)

      do 10 ispin=1,iuhf+1
C
C Watson
C   Grab the AO density and transform to MO basis
C
        call getrec(20,'JOBARC','SCFEVCA0',nbas*nbas*iintfp,scr(i000))
        call getrec(20,'JOBARC','AOOVRLAP',nbas*nbas*iintfp,scr(i010))

        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I010),  nbas,  ! Overlap
     +                           SCR(I000),  nbas,  ! Coeff
     +                     0.0D0,SCR(I020),  nbas ) ! SC

        CALL XGEMM('T','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I000),  nbas,  ! Ct
     +                           SCR(I020),  nbas,  ! SC
     +                     0.0D0,DENS     ,  nbas ) ! Ct SC

C        CALL OUTPUT(DENS,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1

        call getrec(20,'JOBARC',DLABEL,nbas*nbas*iintfp,dens)

        if (dens_diff) then
          call getrec(20,'JOBARC','EXCDEN'//iroot,
     +                       nbas*nbas*iintfp,edens)

          DO i = 1, nbas
          DO j = 1, nbas
             dens (i,j) = edens (i,j) - dens (i,j)
          ENDDO
          ENDDO
        endif

        CALL  DCOPY (NBAS*NBAS,DENS,1,SCR(I000),1)

        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I000),  nbas,
     +                           SCR(I010),  nbas,
     +                     0.0D0,SCR(I020),  nbas )

        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,SCR(I010),  nbas,
     +                           SCR(I020),  nbas,
     +                     0.0D0,SCR(I000),  nbas )

        CALL  DCOPY (NBAS*NBAS,SCR(I000),1,DENS,1)

        call ao2mo2(scr(i000),dens,nlorb,scr(i030),nbas,nbas,ispin)

        DO imo = 1,NBAS
        DO jmo = 1,NBAS
           DVAL = DABS (DENS (imo,jmo))
           IF (DVAL .LT. DENS_THRESH) DENS (imo,jmo) = 0.0D0
        ENDDO
        ENDDO

        WRITE (*,*) ' CSDSC before eig! '
        CALL OUTPUT(DENS,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
C
C
C     Diagonalize MO basis density
C
C
        call zero(nlorb,NBAS*NBAS)
        WRITE (*,*) ' Nat. Orbs. MO x MO ! '
        call eig (dens,nlorb,ijunk,nbas,-1)
        CALL OUTPUT(nlorb,1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1
        call dcopy(nbas,dens,nbas+1,iocc,1)
        nelec = 0.0D0
        DO imo = 1, nbas
           nelec = nelec + iocc(imo)
        ENDDO
        WRITE (*,*) ' Trace of diagonalized MO density matrix - ',nelec

        call getrec(20,'JOBARC','SCFEVCA0',nbas*nbas*iintfp,dens)
        CALL XGEMM('N','N',nbas,nbas,nbas,
     +                     1.0D0,DENS,   nbas,
     +                           NLORB,  nbas,
     +                     0.0D0,SCR(I010),    nbas )

        WRITE (*,*) ' Nat. Orbs. AO x MO ! '
        CALL OUTPUT(SCR(I010),1,NBAS,1,NBAS,NBAS,NBAS,1) ! = 1


        call getrec(20,'JOBARC','CMP2ZMAT',nao*nbas*iintfp,scr(i000))
        WRITE (*,*) ' CMP2ZMAT ! '
        CALL OUTPUT(scr(i000),1,NAO,1,NBAS,NAO,NBAS,1) ! = 1
        call xgemm('n','n',nao,nbas,nbas,one,scr(i000),nao,scr(i010),
     &    nbas,zilch,orb,nao)


        WRITE (*,*) ' Nat. Orbs. NAOBFNS x MO ! '
        CALL OUTPUT(orb,1,NAO,1,NBAS,NAO,NBAS,1) ! = 1

        call getrec(-1,'JOBARC',cscfener,nbas*iintfp,ener)



        do 80 imo=1,nbas
          write(iunit,90)' Ene=   ',ener(imo)
90        format(A,F8.4)
          write(iunit,100)' Spin= ' // sptype(ispin)
100       format(A)
          write(iunit,110)'  Occup=   ',iocc(imo)
110       format(A,F12.6)
          do 120 iao=1,nao
            write(iunit,130)iao,'  ',orb(iao,imo)
130         format(I4,A,F12.6)
120       continue
80      continue 
10    continue
      return
      end




