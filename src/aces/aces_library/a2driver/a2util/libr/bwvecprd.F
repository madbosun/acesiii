
C COMPUTES C(I) = a(i) / (1./b(i) + denomshift)
C FOR FIRST N ELEMENTS OF VECTORS A AND B.

      SUBROUTINE BWVECPRD(denomshift,A,B,C,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N),C(N)
#include "bwcc.com"
#ifdef _ASSERT
      if (n.lt.0) then
         print *, '@BWVECPRD: Assertion failed.'
         print *, '           n = ',n
         call errex
      end if
#endif /* _ASSERT */
cjp b(i) was already 1/(e+e-e-e) and in original routine c was a*b
      DO I = 1, N
c         if (bwgossip) write(*,*) '@BWVECPRD: ',i,b(i),a(i)
         b(i) = 1.0d0 / ( denomshift + (1.0d0/b(i)) )
cjp intruder prevention
         if (abs(b(i)).gt.intruder) then
            c(i) = 0.0d0
         else
            c(i) = a(i)* b(i)
         end if
      END DO
      RETURN
      END
