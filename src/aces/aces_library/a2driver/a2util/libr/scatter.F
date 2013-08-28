
C EMULATOR OF CRAY SCILIB ROUTINE. SCATTERS A VECTOR (B) INTO ANOTHER
C VECTOR (A) ACCORDING TO AN INDEX VECTOR (INDEX).

c The UNICOS man page for scatter is appended to this file.

      subroutine scatter(n,a,index,b)
      integer n,index(n),i
      double precision a(*),b(n)
      do i = 1, n
         a(index(i)) = b(i)
      end do
      return
      end

c NAME
c      SCATTER - Scatters a vector into another vector
c
c SYNOPSIS
c      CALL SCATTER(n,a,index,b)
c
c DESCRIPTION
c      SCATTER is defined as follows:
c
c           a(j(i)) = b(i) where i = 1,...,n
c
c      This routine has the following arguments:
c
c      n      Integer.  (input)
c             Number of elements in arrays index and b (not in a).
c
c      a      Real or integer array of dimension max(index(i): i=1,...,n).
c             (output)
c             Contains the result vector.
c
c      b      Real or integer array of dimension n.  (input)
c             Contains the source vector.
c
c      index  Integer array of dimension n.  (input)
c             Contains the vector of indices.
c
c      The Fortran equivalent of this routine is as follows:
c
c                DO 100 I=1,N
c                   A(INDEX(I))=B(I)
c            100 CONTINUE
c
c CAUTIONS
c      You should not use this routine on systems that have Compress-Index
c      Gather-Scatter (CIGS) hardware, because it will degrade performance.

