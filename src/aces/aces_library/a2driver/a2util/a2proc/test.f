
c This routine grabs a record from JOBARC and compares it to data in test files.
c Recognized datatypes: i/nt/eger, d/ouble, f/loat, r/eal (technically,
c anything not an integer is processed as a double)




      subroutine test(args,numfiles)
      implicit none

c ARGUMENT LIST
      integer numfiles
      character*80 args(numfiles)

c INTERNAL VARIABLES
      integer iFile, iStat
      integer iFirst, iLast
      logical bWarn, bCont
      character*80 sz
      character*80 szRecName
      character*80 szDataType
      character*1 czSpace, czTab

c EXTERNAL FUNCTIONS
      logical leq
      integer fnblnk, linblnk
      external leq, fnblnk, linblnk

c COMMON BLOCKS
      integer            iErrExit
      common /exit_stat/ iErrExit
      save   /exit_stat/


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





c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end

c ----------------------------------------------------------------------

      if (numfiles.lt.1) then
         print *, '@TEST: The test module requires at least 1 argument'
         print '()'
         print *, '  xa2proc test <file> [ <file> [ <file> [ ... ]]]'
         call aces_exit(1)
      end if

      iErrExit = 0

c   o constants
      czSpace = achar(32)
      czTab   = achar(9)

c   o loop over the file arguments
      do iFile = 1, numfiles
         iLast = linblnk(args(iFile))
c         print *, '@TEST: processing file "',args(iFile)(1:iLast),'"'

c      o open the file
         open(unit=10,file=args(iFile)(1:iLast),
     &        form='FORMATTED',status='OLD',err=9999,iostat=iStat)
         rewind(10)
         bWarn = .true.

c      o loop over records and data until EOF or bCont=F
         read(10,'(a)',end=100,err=9999) sz
         bCont = .true.
         do while (bCont)
         if (0.eq.fnblnk(sz)) then
            read(10,'(a)',end=100,err=9999) sz
         else
            bWarn = .false.

c         o read the data type
            iFirst = 1
            do while ((sz(iFirst:iFirst).eq.czSpace.or.
     &                 sz(iFirst:iFirst).eq.czTab).and.
     &                iFirst.le.80)
               iFirst = iFirst + 1
            end do
            iLast = iFirst
            do while ((sz(iLast:iLast).ne.czSpace.and.
     &                 sz(iLast:iLast).ne.czTab).and.
     &                iLast.lt.80)
               iLast = iLast + 1
            end do
            if (iFirst.ne.81) then
               if (sz(iLast:iLast).eq.czSpace.or.
     &             sz(iLast:iLast).eq.czTab) iLast = iLast - 1
               szDataType = sz(iFirst:iLast)
            end if

c         o read the record name
            iFirst = iLast + 1
            do while ((sz(iFirst:iFirst).eq.czSpace.or.
     &                 sz(iFirst:iFirst).eq.czTab).and.
     &                iFirst.le.80)
               iFirst = iFirst + 1
            end do
            iLast = iFirst
            do while ((sz(iLast:iLast).ne.czSpace.and.
     &                 sz(iLast:iLast).ne.czTab).and.
     &                iLast.lt.80)
               iLast = iLast + 1
            end do
            if (iLast.gt.80) then
               print *, '@TEST: string overflow while reading rec name'
               call aces_exit(1)
            end if
            if (sz(iLast:iLast).eq.czSpace.or.
     &          sz(iLast:iLast).eq.czTab) iLast = iLast - 1
            szRecName = sz(iFirst:iLast)

c            print *, 'DEBUG: ',szDataType(1:1),' "',szRecName(1:8),'"'
            if (szDataType(1:1).eq.'i'.or.szDataType(1:1).eq.'I') then
c            o integers
               call test_int(10,sz,bCont,
     &                       szRecName(1:8),icore(i0),icrsiz)
            else
c            o doubles
               call test_dbl(10,sz,bCont,
     &                       szRecName(1:8),icore(i0),icrsiz/iintfp)
            end if

c        end if (sz.ne.' ')
         end if
c        end do while (bCont)
         end do

c      o close the file
 100     close(10,status='KEEP')
         if (bWarn) then
            print *, '@TEST: WARNING: Test file "',args(iFile)(1:iLast),
     &               '" has no test data.'
         end if
      end do

      if (iErrExit.ne.0) call aces_exit(iErrExit)
      print *, '@TEST: All test records pass.'

      return

 9999 call aces_io_error('TEST',10,iStat)

c     end subroutine test
      end

