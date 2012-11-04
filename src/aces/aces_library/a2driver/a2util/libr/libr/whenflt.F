
C EMULATOR OF CRAY SCILIB ROUTINE. LOCATES ALL ELEMENTS OF A DOUBLE
C VECTOR (X) LESS THAN A PARTICULAR VALUE (TARGET).

c The UNICOS man page for whenflt is appended to this file.

      subroutine whenflt(n,x,incx,target,index,nn)
      implicit none
      integer n, incx, index(*), nn
      double precision x(*), target
      integer i, ndx

#ifdef _ASSERT
      if (n.lt.0) then
         print *, '@WHENFLT: Assertion failed.'
         print *, '          n = ',n
         call errex
      end if
#endif /* _ASSERT */

      if (incx.ge.0) then
         if (incx.eq.1) then
            nn = 0
            do i = 1, n
               if (x(i).lt.target) then
                  nn = nn + 1
                  index(nn) = i
               end if
            end do
         else
            nn = 0
            ndx = 1
            do i = 1, n
               if (x(ndx).lt.target) then
                  nn = nn + 1
                  index(nn) = i
               end if
               ndx = ndx + incx
            end do
         end if
c     else if (incx.lt.0) then
      else
         print *, '@WHENFLT: negative increments must be coded'
         call errex
c     end if (incx.ge.0)
      end if

      return
      end

c NAME
c      WHENFLT, WHENFLE, WHENFGT, WHENFGE - Searches a real vector for all
c      elements with a specified logical relationship to a real target
c
c SYNOPSIS
c      CALL WHENFLT(n,x,incx,target,index,nn)
c      CALL WHENFLE(n,x,incx,target,index,nn)
c      CALL WHENFGT(n,x,incx,target,index,nn)
c      CALL WHENFGE(n,x,incx,target,index,nn)
c
c DESCRIPTION
c      WHENFLT searches a real vector for all elements that are less than a
c      real target.
c
c      WHENFLE searches a real vector for all elements that are less than or
c      equal to a real target.
c
c      WHENFGT searches a real vector for all elements that are greater than
c      a real target.
c
c      WHENFGE searches a real vector for all elements that are greater than
c      or equal to a real target.
c
c      When using the Cray Fortran compiler UNICOS systems, all arguments
c      must be of default kind unless documented otherwise.  On UNICOS
c      systems, default kind is KIND=8 for integer, real, complex, and
c      logical arguments.
c
c      These routines have the following arguments:
c
c      n      Integer.  (input)
c             Number of elements of the array to be searched.
c
c      x      Real array of dimension (n-1)* |incx| + 1.  (input)
c             Array x contains the vector to be searched.
c
c      incx   Integer.  (input)
c             Increment between elements of the searched array.
c
c      target Real.  (input)
c             Value for which array is searched.
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
c      for all k = 1, 2, ..., nn.

