
C CALCULATES THE DISTANCE BETWEEN TWO POINTS IN CARTESIAN SPACE

      Double Precision FUNCTION DIST(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),B(3)
cYAU - I warn you, expanding all of these out does not ensure Z is positive.
      Z=  (A(1)-B(1))**2
      Z=Z+(A(2)-B(2))**2
      Z=Z+(A(3)-B(3))**2
      DIST = SQRT(Z)
      RETURN
      END

