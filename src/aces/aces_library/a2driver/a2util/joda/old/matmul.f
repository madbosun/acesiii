      SUBROUTINE MATMUL(A,B,C,NA,NB,NC,NTA,NTB,NTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Maximum number of atoms currently allowed
      PARAMETER (MXATMS = 100)
      DIMENSION B(NA,NB),C(NB,NC),A(NA,NC),SCRATCH(3*MxAtms,3*MxAtms)
      DO 10 I=1,NTA
      DO 10 J=1,NTC
      Z=0.D0
      DO 20 K=1,NTB
20    Z=Z+B(I,K)*C(K,J)
      SCRATCH(I,J)=Z
10    CONTINUE
      DO 15 I=1,NTA
      DO 15 J=1,NTC
15    A(I,J)=SCRATCH(I,J)
      RETURN
      END
