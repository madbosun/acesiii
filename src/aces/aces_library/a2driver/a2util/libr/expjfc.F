      SUBROUTINE EXPJFC(JFC,JFC2,NCOORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION JFC,JFC2
      DIMENSION JFC(NCOORD,NCOORD),JFC2(3*NCOORD,3*NCOORD)

#ifdef _ASSERT
      if (ncoord.lt.0) then
         print *, '@EXPJFC: Assertion failed.'
         print *, '         ncoord = ',ncoord
         call errex
      end if
#endif /* _ASSERT */

      CALL ZERO(JFC2,9*NCOORD*NCOORD)

      INEW = 1
      DO I = 1, NCOORD
         JNEW = 0
         DO J = 1, NCOORD
            INEW1 = INEW
            JNEW  = JNEW + 1
            JFC2(INEW1,JNEW) = JFC(I,J)
            JNEW  = JNEW  + 1
            INEW1 = INEW1 + 1
            JFC2(INEW1,JNEW) = JFC(I,J)
            JNEW  = JNEW  + 1
            INEW1 = INEW1 + 1
            JFC2(INEW1,JNEW) = JFC(I,J)
         END DO
         INEW = INEW + 3
      END DO

      RETURN
      END
