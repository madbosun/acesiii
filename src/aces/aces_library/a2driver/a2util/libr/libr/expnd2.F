
C This routine expands a triangularly packed vector of numbers
C into a square matrix.
C
C  WPACK((NDIM*(NDIM+1))/2) ==> WFULL(NDIM,NDIM)

      SUBROUTINE EXPND2(WPACK,WFULL,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WPACK((NDIM*(NDIM+1))/2),WFULL(NDIM,NDIM)
#ifdef _ASSERT
      if (ndim.lt.0) then
         print *, '@EXPND2: Assertion failed.'
         print *, '         ndim = ',ndim
         call errex
      end if
#endif /* _ASSERT */
      ITHRU = 0
      DO I = 1, NDIM
         DO J = 1, I
            ITHRU = ITHRU + 1
            WFULL(I,J) = WPACK(ITHRU)
            WFULL(J,I) = WPACK(ITHRU)
         END DO
      END DO
      RETURN
      END
