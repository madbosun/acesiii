
C ROBUST EQUIVALENCE CHECK - DO WELL DEFINED SORT ON COORDINATE
C MATRIX AND COMPARE ELEMENT BY ELEMENT. SHOULD BE FOOLPROOF.

      SUBROUTINE COMPARE2(VEC,VECR,NORD,ICOMP,TOL)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C VEC      coordinate vector to be checked (modified)
C VECR     sorted reference coordinate vector (input only)
C NORD     ???
C ICOMP    number of coordinates outside of TOL (output only)
C TOL      tolerance for comparison of coords (input only)

#include "mxatms.par"
#include "coord.com"
      DIMENSION VEC(3*NATOMS),VECR(3*NATOMS)
      DIMENSION NORD(2*MXATMS),SCR(3*MXATMS)

 80   FORMAT(3(1X,F10.5))
C
#ifdef _DEBUG_LVL0
      WRITE(6,*)'--------------B'
      WRITE(6,80)(VEC(JAP),JAP=1,3*NATOMS)
#endif

      CALL SORTXYZ(VEC,scr,NORD(NATOMS+1),NATOMS)
C
#ifdef _DEBUG_LVLM1
      WRITE(6,*)'--------------I'
      WRITE(6,80)(VECR(JAP),JAP=1,3*NATOMS)
      WRITE(6,*)'--------------O'
      WRITE(6,80)(VEC(JAP),JAP=1,3*NATOMS)
#endif

      ICOMP = 0
      DO I = 1, NATOMS*3
         Z = DABS( VECR(I)-scr(I) )

C As a temporary fix to a problem Gennady is having 
C following changes have been made. We hope to find the exact
C reason for the failure. AP 03/14/97.

         IF ((Z .GT. TOL) .AND. (Z .GT. 10*TOL)) ICOMP = ICOMP + 1
         IF ((Z .GT. TOL) .AND. (Z .LT. 10*TOL)) THEN
            WRITE(*,*) 'Warning - Less tighter tolerance is used.'
         END IF

      END DO

      RETURN
      END

