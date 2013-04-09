
C THIS ROUTINE ACCEPTS AN INPUT NSIZE x NSIZE SQUARE MATRIX [VFULL]
C AND CONSTRUCTS THE PACKED LOWER TRIANGULAR MATRIX VPACK [LENGTH:
C (VFULL*VFULL+1)/2

      SUBROUTINE SQUEZ2(VFULL,VPACK,NSIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VFULL(NSIZE,NSIZE),VPACK(NSIZE*(NSIZE+1)/2)
#ifdef _ASSERT
      if (nsize.lt.0) then
         print *, '@SQUEZ2: Assertion failed.'
         print *, '         nsize = ',nsize
         call errex
      end if
#endif /* _ASSERT */
      ITHRU = 0
      DO I = 1, NSIZE
         DO J = 1, I
            ITHRU = ITHRU + 1
            VPACK(ITHRU) = VFULL(I,J)
         END DO
      END DO
      RETURN
      END
