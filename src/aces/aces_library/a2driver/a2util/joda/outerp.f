
C SUBROUTINE COMPUTES OUTER PRODUCT OF VECTOR V WITH ITSELF.  IF
C ICOMP IS SET TO ONE, 1-|V><V| IS RETURNED; OTHERWISE IT IS JUST
C |V><V|.

      SUBROUTINE OUTERP(VV,V,NR,NC,ICOMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(NR,NC),VV(NR,NR)
      DO I = 1, NR
         DO J = 1, NR
            T = 0.D0
            DO K = 1, NC
               T = T + V(I,K) * V(J,K)
            END DO
            VV(I,J) = T
         END DO
      END DO
      IF (ICOMP.EQ.0) RETURN
      Z = -1.D0
      CALL xscal(NR*NR,Z,VV,1)
      DO I = 1, NR
         VV(I,I) = 1.D0 + VV(I,I)
      END DO
      RETURN
      END

