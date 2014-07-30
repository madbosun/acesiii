
c A debugging routine that turns dline(izl(1,i):izl(2,i)) to carets (^)
c repeating for j while making the rest spaces.




      subroutine prep_dline_2(dline,izl,i,j)
      implicit none

c     Maximum string length of terminal lines
      INTEGER LINELEN
      PARAMETER (LINELEN=80)

      character*(linelen) dline
      integer izl(2,7), i, j

      integer ndx, k

      character*1 achar, czTmp, czSpace, czCaret

c ----------------------------------------------------------------------

      czCaret = achar(94)
      czSpace = achar(32)

      do ndx = 1, 7
         if (izl(1,ndx).ne.0) then
            if ((ndx.eq.i).or.(ndx.eq.j)) then
               czTmp = czCaret
            else
               czTmp = czSpace
            end if
            do k = izl(1,ndx), izl(2,ndx)
               dline(k:k) = czTmp
            end do
         end if
      end do

      return
      end

