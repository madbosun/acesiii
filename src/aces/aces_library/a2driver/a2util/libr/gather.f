
C EMULATOR OF CRAY SCILIB ROUTINE. GATHERS A VECTOR (B) INTO ANOTHER
C VECTOR (A) ACCORDING TO AN INDEX VECTOR (INDEX).

c The UNICOS man page for gather is appended to this file.

      subroutine gather(n,a,b,index)
      integer n,index(n),i
      double precision a(n),b(*)
      do i = 1, n
         a(i) = b(index(i))
      end do
      return
      end

c NAME
c      GATHER - Gathers a vector from a source vector
c
c SYNOPSIS
c      CALL GATHER(n,a,b,index)
c
c DESCRIPTION
c      GATHER is defined as follows:
c
c           a(i) = b(j(i)) where i = 1,...,n
c
c      This routine has the following arguments:
c
c      n      Integer.  (input)
c             Number of elements in arrays a and index (not in b).
c
c      a      Real or integer array of dimension n.  (output)
c             Contains the result vector.
c
c      b      Real or integer array of dimension max(index(i): i=1,...,n).
c             (input)
c             Contains the source vector.
c
c      index  Integer array of dimension n.  (input)
c             Contains the vector of indices.
c
c      The Fortran equivalent of this routine is as follows:
c
c                DO 100 I=1,N
c                   A(I)=B(INDEX(I))
c            100 CONTINUE
c
c CAUTIONS
c      You should not use this routine on systems that have Compress-Index
c      Gather-Scatter (CIGS) hardware, because it will degrade performance.

