
C EMULATOR OF CRAY SCILIB ROUTINE. LOCATES ALL ELEMENTS OF AN INTEGER
C VECTOR (IA) EQUAL TO A PARTICULAR VALUE (TARGET).

c The UNICOS man page for wheneq is appended to this file.
c NOTE: Since Cray computers always use 64-bit integers, the vendor's
c       wheneq accepts both integers and floats in the array argument.
c       This version only accepts integer arrays.

      subroutine wheneq(n,ia,inc,target,index,nn)
      implicit none
      integer n, ia(*), inc, target, index(*), nn
      integer i, ndx

#ifdef _ASSERT
      if (n.lt.0) then
         print *, '@WHENEQ: Assertion failed.'
         print *, '         n = ',n
         call errex
      end if
#endif /* _ASSERT */

      if (inc.ge.0) then
         if (inc.eq.1) then
            nn = 0
            do i = 1, n
               if (ia(i).eq.target) then
                  nn = nn + 1
                  index(nn) = i
               end if
            end do
         else
            nn = 0
            ndx = 1
            do i = 1, n
               if (ia(ndx).eq.target) then
                  nn = nn + 1
                  index(nn) = i
               end if
               ndx = ndx + inc
            end do
         end if
c     else if (inc.lt.0) then
      else
         print *, '@WHENEQ: negative increments must be coded'
         call errex
c     end if (inc.ge.0)
      end if

      return
      end

c NAME
c      WHENEQ, WHENNE - Searches a vector for all elements equal or not
c      equal to a target
c
c SYNOPSIS
c      CALL WHENEQ(n,x,incx,target,index,nn)
c      CALL WHENNE(n,x,incx,target,index,nn)
c
c DESCRIPTION
c      WHENEQ searches a real or integer vector for all elements equal to a
c      real or integer target.
c
c      WHENNE searches a real or integer vector for all elements not equal
c      to a real or integer target.
c
c      When using the Cray Fortran compiler UNICOS systems, all arguments
c      must be of default kind unless documented otherwise.  On UNICOS
c      systems, default kind is KIND=8 for integer, real, complex, and
c      logical arguments.
c
c      These routines have the following arguments:
c
c      n      Integer.  (input)
c             Number of elements to be searched.
c
c      x      Real or integer array of dimension (n-1)* |incx| + 1.  (input)
c             Array x containing the vector to be searched.
c
c      incx   Integer.  (input)
c             Increment between elements of the searched array.
c
c      target Real or integer.  (input)
c             Value for which to search in the array.
c
c      index  Integer array of dimension n.  (output)
c             Array index contains the indices in the array elements that
c             match target.
c
c      nn     Integer.  (output)
c             Number of values put in the index array.
c
c NOTES
c      When scanning backward (incx < 0), each routine starts at the end of
c      the vector and moves backward, as follows:
c
c      x(1-incx*(n-1)), x(1-incx*(n-2)), ..., x(1)
c
c      The desired values are at:
c
c      x(1+(index(k)-1)*incx) when incx > 0
c      x(1+(index(k)-n)*incx) when incx < 0
c
c      for all k = 1, 2, ..., nn .

