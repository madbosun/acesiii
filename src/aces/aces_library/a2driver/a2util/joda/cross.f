      SUBROUTINE CROSS(A,B,C,IX)
C
C CALCULATES THE (OPTIONALLY) NORMALIZED VECTOR CROSS PRODUCT C=A x B
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),B(3),C(3)
      C(3)=A(1)*B(2)-B(1)*A(2)
      C(2)=-A(1)*B(3)+A(3)*B(1)
      C(1)=A(2)*B(3)-A(3)*B(2)
C
C This block is added to avoid calling NORMAL when
C the vector is null. The null vector is reset to
C vector A. Ajith Perera, 10/2010
C
      D2 = 0.0D0
      Do I = 1, 3
         D2 = D2 + C(i)**2
      Enddo
      D2sqrt = Dsqrt(D2)
      If (D2sqrt .LT. 1.0D-14) Then
      Return
      Else
          IF(IX.EQ.1)CALL NORMAL(C,3)
      Endif
C
      RETURN
      END
