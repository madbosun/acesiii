      DOUBLE PRECISION FUNCTION ZETA(I,J,K)
      IF (I.EQ.J) THEN
         ZETA = 1.0D0
      ELSE IF (I.EQ.K) THEN
         ZETA = -1.0D0
      ELSE
         ZETA = 0.0D0
      END IF
      END
