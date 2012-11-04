
C SUMS LIST (ISPIN,LIST1) AND (ISPIN,LIST2) AND THEN OVERWRITES
C (ISPIN,LIST2) WITH THE SUM.

      SUBROUTINE SUMLST(ISPIN,LIST1,LIST2,Z,T2,NDSSIZ,NDIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T2(NDSSIZ),Z(NDSSIZ,NDIS)
      PARAMETER (ONE=1.0D0)
#ifdef _ASSERT
      if (ndis.lt.0) then
         print *, '@SUMLST: Assertion failed.'
         print *, '         ndis = ',ndis
         call errex
      end if
#endif /* _ASSERT */
      CALL GETLST(Z,1,NDIS,1,ISPIN,LIST1)
      DO I = 1, NDIS
         CALL GETLST(T2,I,1,1,ISPIN,LIST2)
         CALL SAXPY(NDSSIZ,ONE,Z(1,I),1,T2,1)
         CALL PUTLST(T2,I,1,1,ISPIN,LIST2)
      END DO
      RETURN
      END
