      SUBROUTINE TIEDIHED(STEP, NOPT, ITIE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION STEP(NOPT),ITIE(NOPT)
C
      CALL GETREC(20,'JOBARC','TIEDCORD',NOPT, ITIE)

      DO  I=1,NOPT
          IF (ITIE(I).NE.0)THEN
             DO J=I+1,NOPT
                IF (ITIE(J).EQ.ITIE(I))THEN
                    STEP(J)=-STEP(I)
      Write(6, "(1x,a,I4,a,a,F13.8)") "Optimized coordinate ",J, 
     &                             " automatically ","set to ",
     &                              -step(i)
 
                ENDIF 
             ENDDO
         ENDIF
      ENDDO
C
      RETURN
      END 
