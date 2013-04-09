
cSG/JDW 6/6/95
c Change of name for compatibility with newer IBM operating systems (XLF 3.2).

      character*(*) function trmblk(sz)
      implicit none

      character*(*) sz
      integer i, fnblnk, linblnk, iStart, iEnd

      trmblk=' '

      iStart = fnblnk(sz)
      if (iStart.ne.0) then
         iEnd = min(linblnk(sz),iStart-1+len(trmblk))
         trmblk = sz(iStart:iEnd)
      end if

      return
      end

