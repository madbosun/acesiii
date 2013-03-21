
c This routine loads a block of lists into icore-addressable memory.

      subroutine aces_auxcache_ldblk(i0,iCrSiz,
     &                               iLeft0, iLeft1, iLeftInc,
     &                               iRight0,iRight1,iRightInc)
      implicit none

c ARGUMENTS
      integer i0, iCrSiz
      integer iLeft0, iLeft1, iLeftInc, iRight0, iRight1, iRightInc

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer iFam, iGrp, nRows, nCols, nInts, i

c COMMON BLOCKS


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



c auxcache.com : begin

c The auxiliary cache is a programmer-controlled list cache. The dimensions are
c the same as those of the MOIO arrays. The programmer may load lists into icore
c memory and set quikget(?,?) to their icore addresses. When getlst and putlst
c operate on ANY list, quikget is checked to see if the list lives in icore.
c If so, the operation is performed on the in-core data instead of hitting the
c storage file(s).

c WARNING
c    There is no automatic updating of storage files with data in the auxiliary
c cache. If the memory-resident data is altered and must be stored to disk, then
c the quikget value must be destroyed (zeroed) and the actual location of the
c data must then be passed to putlst. The aces_auxcache_flush routine
c systematically stores the quikget values, calls putlst, and restores the
c quikget values.





      integer           quikget(10,500)
      common /auxcache/ quikget
      save   /auxcache/

c auxcache.com : end

c ----------------------------------------------------------------------

      i = 0
c   o assert icore info is sane
      if (i0.eq.0.or.iCrSiz.lt.0) then
         print *, '@ACES_AUXCACHE_LDBLK: Assertion failed.'
         print *, '   i0     = ',i0
         print *, '   iCrSiz = ',iCrSiz
         i = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft0 .lt.1).or.(10.lt.iLeft0 ).or.
     &    (iLeft1 .lt.1).or.(10.lt.iLeft1 ).or.
     &    (iRight0.lt.1).or.(500.lt.iRight0).or.
     &    (iRight1.lt.1).or.(500.lt.iRight1).or.
     &    (iLeftInc.eq.0).or.(iRightInc.eq.0)            ) then
         print *, '@ACES_AUXCACHE_LDBLK: Assertion failed.'
         print *, '   do iFam = ',iRight0,',',iRight1,',',iRightInc
         print *, '   do iGrp = ',iLeft0 ,',',iLeft1 ,',',iLeftInc
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)

c ----------------------------------------------------------------------

c   o load lists
      do iFam = iRight0, iRight1, iRightInc
      do iGrp = iLeft0,  iLeft1,  iLeftInc
         if (quikget(iGrp,iFam).eq.0) then
            nCols = aces_list_cols(iGrp,iFam)
            nRows = aces_list_rows(iGrp,iFam)
            nInts = iintfp*nRows*nCols
            if (0.lt.nInts.and.nInts.lt.iCrSiz) then
c               print *, '@ACES_AUXCACHE_LDBLK: Loading ',iGrp,iFam
               call getlst(icore(i0),1,nCols,1,iGrp,iFam)
               quikget(iGrp,iFam) = i0
               i0     = i0     + nInts
               iCrSiz = iCrSiz - nInts
            end if
         end if
      end do
      end do

      return
c     end subroutine aces_auxcache_ldblk
      end

