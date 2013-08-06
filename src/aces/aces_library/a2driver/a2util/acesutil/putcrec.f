
c This routine writes a character string to the job archive.

c WARNING!
c    The "record length" that gets associated with the logical record is the
c number of whole integers that contain the string. Therefore, calling getrec
c with the query flag may not return the actual length of the string.

c INPUT
c int      XFLAG     : (same as putrec)
c char*(*) SZARCHIVE : (same as putrec)
c char*(*) SZRECNAME : (same as putrec)
c int      ILENGTH   : the substring-length of the record to store
c char*(*) SZSRC     : (same as putrec)









      subroutine putcrec(xFlag,szArchive,szRecName,iLength,szSrc)
      implicit none

c ARGUMENTS
      integer xFlag, iLength
      character*(*) szArchive, szRecName, szSrc(*)

c EXTERNAL FUNCTIONS
      integer iachar

c INTERNAL VARIABLES
      integer i, iNdx, nLeft, iBuf(1024), iRecLen

c ----------------------------------------------------------------------


      i = 0
c   o assert substring fits in iBuf
      if (iLength.gt.(1024*8)) then
         print *, '@PUTCREC: Assertion failed.'
         print *, '   iLength = ',iLength
         print *, '   maximum = ',1024*8
         i = 1
      end if
cYAU: Since szSrc is now an array pointer, we have no idea if iLength is good.
cc   o assert substring is within szSrc
c      if (len(szSrc).lt.iLength) then
c         print *, '@PUTCREC: Assertion failed.'
c         print *, '   iLength = ',iLength
c         print *, '   len(sz) = ',len(szSrc)
c         i = 1
c      end if
c   o assert ichar and iachar are the same
      if (ichar(' ').ne.iachar(' ')) then
         print *, '@PUTCREC: Assertion failed.'
         print *, '   iachar(szSpace) = ',iachar(' ')
         print *, '    ichar(szSpace) = ',ichar(' ')
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)


c ----------------------------------------------------------------------

c   o iRecLen is the number of whole integers that contain szSrc(1:iLength)
      iRecLen = (iLength+8-1)/8

c   o flush the last integer with spaces
      if (iand(iLength,8-1).ne.0) then
         i = ichar(' ')
         call c_memset(iBuf(iRecLen),i,8)
      end if

c   o copy the string to the integer buffer
      call c_memmove(iBuf,szSrc,iLength)

c   o write to the job archive
      call putrec(xFlag,szArchive,szRecName,iRecLen,iBuf)

      return
c     end subroutine putcrec
      end

