
C THIS ROUTINE EXPANDS A COMPRESSED ARRAY TO A FULL ANTISYMMETRIC ARRAY

      SUBROUTINE EXPND3(WPACK,WFULL,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      DIMENSION WPACK((NDIM*(NDIM-1))/2),WFULL(NDIM,NDIM)
      DIMENSION WPACK(*),WFULL(NDIM,NDIM)
#ifdef _ASSERT
      if (ndim.lt.0) then
         print *, '@EXPND3: Assertion failed.'
         print *, '         ndim = ',ndim
         call errex
      end if
#endif /* _ASSERT */
      ITHRU = 0
      DO I = 2, NDIM
         DO J = 1, I-1
            ITHRU = ITHRU + 1
            WFULL(I,J) = -WPACK(ITHRU)
            WFULL(J,I) =  WPACK(ITHRU)
         END DO
      END DO
      DO I = 1, NDIM
         WFULL(I,I) = 0.0d0
      END DO
      RETURN
      END
