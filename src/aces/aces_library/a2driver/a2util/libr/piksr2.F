
C SIMPLE SORTER FROM NUMERICAL RECIPES.

      SUBROUTINE PIKSR2(N,ARR,NLIST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARR(N),NLIST(N)
      if (n.lt.2) return
      DO J = 2, N
         A = ARR(J)
         NLST = NLIST(J)
         DO I = J-1, 1, -1
            IF (ARR(I).LE.A) GOTO 10
            ARR(I+1) = ARR(I)
            NLIST(I+1) = NLIST(I)
         END DO
         I = 0
   10    ARR(I+1) = A
         NLIST(I+1) = NLST
      END DO
      RETURN
      END
