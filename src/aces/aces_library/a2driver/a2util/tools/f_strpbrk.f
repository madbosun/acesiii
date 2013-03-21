
c This routine returns the first index in sz of any single char appearing in
c charset. WARNING: Since Fortran pads the end of all strings with spaces,
c be sure to pass only the relevant substring.




      integer function f_strpbrk(sz,charset)
      implicit none

      character*(*) sz, charset

      character*(256)      sz2
      character*(32) charset2

      character*1 achar
      integer f_strpbrk_core

      if (len(sz).ge.256) then
         print *, '@F_STRPBRK: The sz buffer is too small ',
     &            'to contain the input string.'
         print *, '                  Recompile with at least ',
     &            len(sz)+1, ' characters in the buffer.'
         print *, '                  (Currently ',256,
     &            ' characters.)'
         call c_exit(1)
      end if
      if (len(charset).ge.32) then
         print *, '@F_STRPBRK: The charset buffer is too small ',
     &            'to contain the input string.'
         print *, '                  Recompile with at least ',
     &            len(charset)+1, ' characters in the buffer.'
         print *, '                  (Currently ',32,
     &                                                   ' characters.)'
         call c_exit(1)
      end if

      sz2      = sz//achar(0)
      charset2 = charset//achar(0)

c      print *, 'sz is ',len(sz),' chars long'
c      print *, 'charset is ',len(charset),' chars long'

      f_strpbrk = f_strpbrk_core(sz2,charset2)

      return
      end

