
C CALCULATES THE VECTOR V BETWEEN CARTESIAN POINTS A AND B

      SUBROUTINE VEC(A,B,V,IX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),B(3),V(3)
      V(1)=B(1)-A(1)
      V(2)=B(2)-A(2)
      V(3)=B(3)-A(3)
C
C This block is added to avoid calling NORMAL when
C the vector is null, Ajith Perera, 10/2010
C
      D2 = 0.0D0
      Do I = 1, 3
         D2 = D2 + V(i)**2
      Enddo
      D2sqrt = Dsqrt(D2)
      If (D2sqrt .LT. 1.0D-14) Then
         Write(6,"(a)") "@-VEC: The null vector is returned"
         If (Ix .NE.0) Call Errex
      Endif
C
      IF(IX.EQ.1)CALL NORMAL(V,3)
      RETURN
      END 
