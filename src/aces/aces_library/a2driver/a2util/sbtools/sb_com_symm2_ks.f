
      subroutine sb_com_symm2_ks
      implicit none



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











      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/




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



      double precision nucrep
      integer length,ijunk,irrep,sqroff,trioff,n,tri,sqr,i,orb,iunit
      character*80 iiiifil
      character*8 title(24)
      logical bExist
      call gfname('IIII',iiiifil,length)
      inquire(file=iiiifil,exist=bExist)
      if (.not.bExist) then
         write(*,*) '@SB_COM_SYMM2_KS: No IIII integral file.'
         call errex
      end if
      iunit = 10
      open(iunit,file=iiiifil,form='unformatted',access='sequential')
      read(iunit) title,nirrep
      rewind(iunit)
      read(iunit) title,nirrep,(numbasir(i),i=1,nirrep),nucrep
      close(iunit,status='keep')
      call getrec(0,'JOBARC','NUCREP  ',length,n)
      if (length.eq.0) call putrec(1,'JOBARC','NUCREP  ',iintfp,nucrep)
      call putrec(1,'JOBARC','NUMBASIR',nirrep,numbasir)
      sqroff=1
      trioff=1
      maxirrsqr=0
      maxirrtri=0
      orb=1
      do irrep=1,nirrep

         n=numbasir(irrep)
         irrorboff(irrep)=orb
         orb=orb+n

         sqr=n*n
         tri=n*(n+1)/2

         irrsqrlen(irrep)=sqr
         irrsqroff(irrep)=sqroff
         sqroff=sqroff+sqr
         maxirrsqr=max(maxirrsqr,sqr)

         irrtrilen(irrep)=tri
         irrtrioff(irrep)=trioff
         trioff=trioff+tri
         maxirrtri=max(maxirrtri,tri)

      end do
      irrorboff(nirrep+1)=orb
      irrsqrtot=sqroff-1
      irrtritot=trioff-1

      call getrec(-1,'JOBARC','OCCUPYA0',nirrep,occup(1,1))
      call getrec(-1,'JOBARC','OCCUPYB0',nirrep,occup(1,2))
      totocc(1)=0
      totocc(2)=0
      do i=1,nirrep
         totocc(1)=totocc(1)+occup(i,1)
         totocc(2)=totocc(2)+occup(i,2)
      end do
      if (iuhf.ne.1) totocc(2)=totocc(1)
      totocca=totocc(1)
      totoccb=totocc(2)
      return
c     end subroutine sb_com_symm2_ks
      end

