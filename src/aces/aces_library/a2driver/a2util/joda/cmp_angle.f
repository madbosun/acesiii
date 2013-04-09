      SUBROUTINE CMP_ANGLE(ILFTEND, IMIDLE, IRHTEND, CARTCOORD,
     &                     VLUANGLE, NRATMS) 
C
C This is a simple routine to calculate the angle between two vectors,
C ILFTEND and IRHTEND and the center is IMIDLE. The angle is retrun
C in VLUANGLE.
C      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION CARTCOORD(3*NRATMS), VEC1(3), VEC2(3)


      CALL VEC(CARTCOORD(IMIDLE*3 - 2), CARTCOORD(ILFTEND*3 - 2), 
     &         VEC1, 1)
C
      CALL VEC(CARTCOORD(IMIDLE*3 - 2), CARTCOORD(IRHTEND*3 - 2), 
     &         VEC2, 1)

      VLUANGLE = ANGLE(VEC1, VEC2, 3)
C
      RETURN
      END
