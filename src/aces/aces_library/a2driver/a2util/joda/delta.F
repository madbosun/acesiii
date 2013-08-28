      FUNCTION DELTA(I,J)
      DOUBLE PRECISION DELTA
      INTEGER I,J
      IF (I.EQ.J) THEN
         DELTA = 1.0D0
      ELSE
         DELTA = 0.0D0
      END IF
      RETURN
      END
