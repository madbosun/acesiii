
C THIS ROUTINE CONVERTS FROM THE MBPT/CC DENOMINATOR
C
C     1/(Fii + Fjj - Faa - Fbb) 
C
C TO THE CI DENOMINATOR
C
C     1/(Fii + Fjj - Faa - Fbb + Ecorr)

      SUBROUTINE CIDENOM(N,ECORR,VEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(1)
c old
c      DO 10 I=1,N
c       X1=VEC(I)
c       X2=1./X1
c       X2=X2+ECORR
c       VEC(I)=1.0/X2
c10    CONTINUE
c new
      if (n.lt.1) return
      DO I = 1, N
         X = 1.0d0 + ( Ecorr * VEC(I) )
         X = 1.0d0 / X
         VEC(I) = VEC(I) * X
      END DO
c end
      RETURN
      END
