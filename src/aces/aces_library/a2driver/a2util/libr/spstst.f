
C RETURNS THE SPARSITY LEVEL OF THE FIRST LEN ELEMENTS OF VEC, EXPRESSED
C AS A DECIMAL FRACTION. TOLERANCE IS SET BY TOL.

      SUBROUTINE SPSTST(VEC,LEN,TOL,SPARSE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(LEN)

      if (len.lt.0) then
         print *, '@SPSTST: Assertion failed.'
         print *, '         len = ',len
         call errex
      end if

      IZRO = 0
      DO I = 1, LEN
         IF (DABS(VEC(I)).LT.TOL) IZRO = IZRO + 1
      END DO
      SPARSE = FLOAT(IZRO) / FLOAT(LEN)
      RETURN
      END
