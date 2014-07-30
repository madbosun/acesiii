
C THIS ROUTINE ACCEPTS VECTORS RV AND IV AND PERFORMS THE
C FOLLOWING OPERATION: IT LEAVES ENTRIES 1-ISKIP UNCHANGED,
C AND OVERWRITES RV(J) AND IV(J) WITH RV(J-1) AND IV(J-1)
C RESPECTIVELY FOR ALL J-1 .GE. ISKIP.  THE VECTORS BOTH HAVE
C LENGTH NUMEL.  RV IS DOUBLE PRECISION AND IV IS INTEGER.
C THIS ROUTINE IS A DEPENDENT OF SCANVC.

      SUBROUTINE SHIFTV(RV,IV,ISKIP,NUMEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RV(NUMEL),IV(NUMEL)
#ifdef _ASSERT
      if (iskip.lt.1) then
         print *, '@SHIFTV: Assertion failed.'
         print *, '         iskip = ',iskip
         call errex
      end if
      if (iskip.gt.numel) then
         print *, '@SHIFTV: Assertion failed.'
         print *, '         numel = ',numel
         print *, '         iskip = ',iskip
         call errex
      end if
#endif /* _ASSERT */
      if ((numel.gt.1).and.(iskip.lt.numel)) then
         DO I = NUMEL-1, ISKIP, -1
            RV(I+1) = RV(I)
            IV(I+1) = IV(I)
         END DO
      end if
      RETURN
      END
