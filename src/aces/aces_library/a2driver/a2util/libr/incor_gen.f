
C This routine loads lists into ICORE(I0), increments I0, and decrements
C ICRSIZ. I0 and ICRSIZ do not have to point to the /ISTART/ common
C block, but they must be relative to ICORE(1) in the blank common
C block. The logic is also reentrant, so it will not load lists that are
C already in memory.




































      SUBROUTINE INCOR_GEN(I0,ICRSIZ,IUHF)
      IMPLICIT NONE
      INTEGER I0,ICRSIZ,IUHF


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
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
      integer           quikget(10,500)
      common /auxcache/ quikget
      integer iLim, iGrp, iFam, nRows, nCols, nDbls, nInts
      integer  aces_list_rows, aces_list_cols
      external aces_list_rows, aces_list_cols

      IF (IFLAGS(35).EQ.0) RETURN

      IF (IFLAGS(35).EQ.6) THEN
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,1,500,1)
         RETURN
      END IF

      iLim = 0
      IF (IFLAGS(35).EQ.1) THEN
c      o NOABCD
         iLim=max(iLim,irpdpd(1,13)*irpdpd(1,11))
         iLim=max(iLim,irpdpd(1,13)*irpdpd(1,18))
         iLim=max(iLim,irpdpd(1,19)*irpdpd(1,9))
         iLim=max(iLim,irpdpd(1,20)*irpdpd(1,10))
      ELSE IF (IFLAGS(35).EQ.4) THEN
c      o NOABCI
         iLim=max(iLim,irpdpd(1,9)*irpdpd(1,10))
         iLim=max(iLim,irpdpd(1,9)*irpdpd(1,9))
         iLim=max(iLim,irpdpd(1,10)*irpdpd(1,10))
      ELSE IF (IFLAGS(35).EQ.5) THEN
c      o NOABIJ
         iLim=max(iLim,irpdpd(1,14)*irpdpd(1,18))
         iLim=max(iLim,irpdpd(1,14)*irpdpd(1,11))
         iLim=max(iLim,irpdpd(1,21)*irpdpd(1,16))
         iLim=max(iLim,irpdpd(1,22)*irpdpd(1,17))
      END IF
      IF (iLim.eq.0) RETURN
      iLim = iintfp*iLim

      do iFam = 1, 500
      do iGrp = 1, 10
         if (quikget(iGrp,iFam).eq.0) then
            nCols = aces_list_cols(iGrp,iFam)
            nRows = aces_list_rows(iGrp,iFam)
            nInts = iintfp*nRows*nCols
            if (0.lt.nInts.and.nInts.le.iLim.and.nInts.lt.iCrSiz) then
c               print *, '@INCOR_GEN: Loading ',iGrp,iFam
               call getlst(icore(i0),1,nCols,1,iGrp,iFam)
               quikget(iGrp,iFam) = i0
               i0     = i0     + nInts
               iCrSiz = iCrSiz - nInts
            end if
         end if
      end do
      end do

      RETURN
      END

