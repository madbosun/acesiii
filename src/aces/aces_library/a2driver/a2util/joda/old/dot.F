      DOUBLE PRECISION FUNCTION DOT(A,B,N)
      INTEGER N, I
      DOUBLE PRECISION A(N),B(N)
      DOT=0.0D0
      DO I=1,N
         DOT=DOT+A(I)*B(I)
      END DO
      RETURN
      END
