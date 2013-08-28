      DOUBLE PRECISION FUNCTION TRACE(A,N)
      IMPLICIT NONE
      INTEGER N, INC, I, ADD
      DOUBLE PRECISION A(N*N)
      TRACE = A(1)
      INC = N + 1
      ADD = 1
      DO I = 1, N-1
         ADD = ADD + INC
         TRACE = TRACE + A(ADD)
      END DO
      RETURN
      END 
