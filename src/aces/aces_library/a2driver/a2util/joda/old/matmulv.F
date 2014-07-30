      SUBROUTINE MATMULV(A,B,C,NA,NB,NC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NC,NA),B(NB,NA),C(NB,NC)
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      CALL XGEMM('T','N',NC,NA,NB,ONE,C,NB,B,NB,ZILCH,A,NC)
c      DO 10 I=1,NA
c      DO 10 K=1,NC
c      Z=0.D0
c      DO 20 J=1,NB
c      Z=Z+B(J,I)*C(J,K)
c20    CONTINUE
c      A(K,I)=Z
c10    CONTINUE
       RETURN 
      END
