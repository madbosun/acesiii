
      SUBROUTINE COMPKMAT(DERBMAT,NRATMS,TOTREDNCO,FI,HC,
     &                    DIFTEMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER TOTREDNCO
      DOUBLE PRECISION KMAT

#include "mxatms.par"

C KMAT is used as temporary storge while building DIFTEMP and it is of
C size 9*MXATMS*MXATMS (so it needs to be managed dynamically).

      DIMENSION DERBMAT(3*NRATMS,3*NRATMS*TOTREDNCO),
     &          FI(TOTREDNCO),HC(3*NRATMS,3*NRATMS),
     &          KMAT(3*MXATMS,3*MXATMS),
     &          DIFTEMP(3*NRATMS,3*NRATMS)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      CALL ZERO(KMAT,9*NRATMS*NRATMS)

      LENDBMAT = 9*NRATMS*NRATMS*TOTREDNCO
      CALL GETREC(20,'JOBARC','DERBMAT',IINTFP*LENDBMAT,DERBMAT)

      DO I=1, TOTREDNCO
         IND=3*NRATMS*(I-1)+1
C
#ifdef _DEBUG_LVL0
         Write(6, "(a,4F10.5)"),'The gradients Compkmat: ',
     &        (FI(J), J=1, TOTREDNCO)
#endif
         CALL DSCAL(9*NRATMS*NRATMS,FI(I),DERBMAT(1,IND),1)
      END DO

      DO I=1,3*NRATMS
         DO J=1,3*NRATMS
            DO K=1,TOTREDNCO
               JTMP = (K-1)*3*NRATMS+J
               KMAT(I,J) = KMAT(I,J)+ DERBMAT(I,JTMP)
            END DO
         END DO
      END DO
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)"),'The  Cartesian Hess'
      CALL OUTPUT(HC,1,3*NRATMS,1,3*NRATMS,3*NRATMS,
     &            3*NRATMS,1)
#endif
C
      DO I=1,3*NRATMS
         DO J=1,3*NRATMS
            DIFTEMP(I,J)=HC(I,J)-KMAT(I,J)
         END DO
      END DO
C
#ifdef _DEBUG_LVL0
      Write(6,*),'In COMPKMAT: The (HC-K) Matrix'
      CALL OUTPUT(DIFTEMP,1,3*NRATMS,1,3*NRATMS,3*NRATMS,
     &            3*NRATMS,1)
#endif

      RETURN
      END

