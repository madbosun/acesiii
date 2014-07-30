      SUBROUTINE CNVRT2RAD(R,NX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      DIMENSION R(NX)
      ATOR=DACOS(-1.D0)/180.D0
      ATOB=0.529177249D0
      DO IX = 8, NX-1, 3
         IF (IX.NE.8) R(IX+1) = R(IX+1)*ATOR
         R(IX) = R(IX)*ATOR
      END DO
      IF (IFLAGS(78).EQ.0) THEN
         DO IX = 4, NX-2, 3
            R(IX) = R(IX)/ATOB
         END DO
      END IF
      RETURN
      END
