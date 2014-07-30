
C ROUTINE ACCEPTS A SYMMETRIC MATRIX V, AND THEN RETURNS THE
C RESULT IN PACKED I.LE.J TRIANGULAR FORM.  THE DIAGONAL IS INCLUDED IF
C IDIAG IS SET TO 0, AND NOT INCLUDED IF IDIAG = 1.

      SUBROUTINE SQUEEZ(V,N,IDIAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N*N)
#ifdef _ASSERT
      if ((idiag.lt.0).or.(idiag.gt.1)) then
         print *, '@SQUEEZ: Assertion failed.'
         print *, '         idiag = ',idiag
         call errex
      end if
#endif /* _ASSERT */
      if (n.lt.1) return
      IX = 0
      DO J = 1+IDIAG, N
         DO I = 1, J-IDIAG
            IX = IX + 1
            INDX = I + (J-1)*N
            V(IX) = V(INDX)
         END DO
      END DO
      RETURN
      END
