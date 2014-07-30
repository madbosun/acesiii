
c This routine returns the index of the smallest number in T(*).
c BEWARE - IMINPOS is 1 if n=0

      INTEGER FUNCTION IMINPOS(N,T,INC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION T(*)
#ifdef _ASSERT
      if (n.lt.0) then
         print *, '@IMINPOS: Assertion failed.'
         print *, '          n = ',n
         call errex
      end if
#endif /* _ASSERT */
      itmp = 1
      XMIN = T(1)
      ndx  = 1 + inc
      DO I = 2, n
         IF (T(ndx).LT.XMIN) THEN
            itmp = i
            XMIN = T(ndx)
         END IF
         ndx = ndx + inc
      END DO
      IMINPOS = itmp
      RETURN
      END
