
c ZEROS OUT THE FIRST LEN ELEMENTS OF INTEGER VECTOR IVEC.

      SUBROUTINE IZERO(IVEC,LEN)
      DIMENSION IVEC(LEN)
      if (len.lt.1) return
      DO I = 1, LEN
         IVEC(I) = 0
      END DO
      RETURN
      END
