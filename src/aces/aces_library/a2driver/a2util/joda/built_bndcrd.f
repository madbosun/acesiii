      SUBROUTINE BULT_BNDCRD(CARTCOORD, BMATRX, DISTAB, ICON1,  
     &                       ICON2, IBNDS, TOTREDNCO, NRATMS)
C         
C Setup the bond stretching B-matrix elements. 
C B(*,i,j)(x,y,z) = (x_i,y_i,z_i - x_j,y_j,z_j)/|Rij|
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER TOTREDNCO
      DIMENSION CARTCOORD(3*NRATMS), BMATRX(TOTREDNCO, 3*NRATMS)
C     
      DISTAB = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
C
      DO 10 IXYZ = 1, 3
C
         DIFFAB = CARTCOORD((3*ICON1 - 3) + IXYZ)
     &             - CARTCOORD((3*ICON2 -3) + IXYZ)
         BMATRX(IBNDS, (3*ICON1 - 3) + IXYZ) =  DIFFAB/DISTAB
         BMATRX(IBNDS, (3*ICON2 - 3) + IXYZ) = -DIFFAB/DISTAB
C
 10   CONTINUE

      RETURN
      END
