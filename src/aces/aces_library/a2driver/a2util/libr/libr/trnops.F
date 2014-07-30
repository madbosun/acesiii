
C TRANSFORMS SYMMETRY OPERATIONS TO AN ORIENTATION SPECIFIED BY
C MATRIX "GOOFY".

      SUBROUTINE TRNOPS(OPS,GOOFY,IORDGP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION OPS(9*IORDGP),GOOFY(9),SCR(9)

#ifdef _ASSERT
      if (iordgp.lt.0) then
         print *, '@TRNOPS: Assertion failed.'
         print *, '         iordgp = ',iordgp
         call errex
      end if
#endif /* _ASSERT */

      DO I = 1, IORDGP
         IADD = -8 + (9*I)
         CALL XGEMM('N','N',3,3,3,
     &              1.0d0,GOOFY,    3,
     &                    OPS(IADD),3,
     &              0.0d0,SCR,      3)
         CALL XGEMM('N','T',3,3,3,
     &              1.0d0,SCR,      3,
     &                    GOOFY,    3,
     &              0.0d0,OPS(IADD),3)
      END DO

      RETURN
      END 

