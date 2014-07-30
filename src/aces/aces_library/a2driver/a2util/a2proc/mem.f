
c This routine changes how much memory ACESII tries to allocate.













































































































































































































































































































































































































































































































































      subroutine mem(args,nargs)
      implicit none

c ARGUMENT LIST
      integer nargs
      character*80 args(nargs)

c INTERNAL VARIABLES
      double precision dMem
      integer i, iFact
      logical bBytes

c COMMON BLOCKS


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

c ----------------------------------------------------------------------

      if (nargs.gt.0) then

c   o find the last non-white character
      i = 80
      do while (args(1)(i:i).eq.' ')
         i = i - 1
      end do

c   o get base unit
      if ((iachar(args(1)(i:i)).lt.iachar('0').or.
     &     iachar(args(1)(i:i)).gt.iachar('9')    ).and.
     &    args(1)(i:i).ne.'.') then
         bBytes=(args(1)(i:i).eq.'B'.or.args(1)(i:i).eq.'b')
         i = i - 1
      end if

c   o get order of magnitude
      iFact=1
      if ((iachar(args(1)(i:i)).lt.iachar('0').or.
     &     iachar(args(1)(i:i)).gt.iachar('9')    ).and.
     &    args(1)(i:i).ne.'.') then
         if (args(1)(i:i).eq.'K'.or.args(1)(i:i).eq.'k') iFact=1024
         if (args(1)(i:i).eq.'M'.or.args(1)(i:i).eq.'m') iFact=1024**2
         if (args(1)(i:i).eq.'G'.or.args(1)(i:i).eq.'g') iFact=1024**3
         i = i - 1
      end if

c   o read the amount
      read(args(1)(1:i),*) dMem

      dMem = dMem*iFact
      if (bBytes) dMem = dMem/iIntLn
      iflags(36) = dMem
      call putrec(1,'JOBARC','IFLAGS',36,iflags)

c     end if (nargs.gt.0)
      end if

      print '()'
      print *, '@MEM: The current value of MEM is'
      i = iflags(36)
      print *, '      ',i,' Integer Words'
      i = i*iIntLn
      print *, '      ',i,' Bytes'
      print *, '      ',i*1.0/1024.0/1024.0,' MegaBytes'
      print *, '      ',i*1.0/1024.0/1024.0/1024.0,' GigaBytes'
      print '()'

      return
c     end subroutine mem
      end

