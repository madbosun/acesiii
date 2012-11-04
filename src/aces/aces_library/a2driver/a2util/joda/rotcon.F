      SUBROUTINE ROTCON(IT,IPRNT,IERR,NOSILENT)
C
C COMPUTES ROTATIONAL CONSTANTS (IN CM-1), OPTIONALLLY
C  PRINTS THEM AND RETURNS AN ERROR MESSAGE IF INERTIA
C  TENSOR PASSED IS NOT DIAGONAL.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION IT(3,3),RC(3),NORD(3)
      LOGICAL NOSILENT
#include "machsp.com"
      DATA FACTOR /60.198D0/
      Z=0.D0
      IBOT=0
      IERR=0
      CALL ZERO(RC,3)
      DO 10 I=1,3
      DO 10 J=I+1,3
       Z=DABS(IT(I,J))+DABS(IT(J,I))+Z
 10   CONTINUE
      IF(Z.GT.1.D-10)THEN
       WRITE(6,100)
 100   FORMAT(T3,' ***PROGRAM ERROR***, Inertia tensor not '
     &,'diagonal in ROTCON.')
       IERR=1
      ELSE
       IBOT=1
      ENDIF
      DO 20 I=1,3
       IF(DABS(IT(I,I)).GT.1.D-12)THEN
        RC(I)=FACTOR/IT(I,I)
       ELSE
        IBOT=IBOT+1
       ENDIF
 20   CONTINUE
      CALL PIKSR2(3,RC,NORD)
C
C Rotational constants are useful things to have in JOBARC.
C
      CALL PUTREC(20, 'JOBARC', 'ROTCONST', 3*IINTFP, RC) 
      IF(IPRNT.NE.0 .AND. NOSILENT)THEN
       WRITE(6,200)
       WRITE(6,300)(RC(J),J=IBOT,3)
 200   FORMAT(T3,' Rotational constants (in cm-1): ')
 300   FORMAT((T3,3(F10.5,5X)))
      ENDIF
      RETURN
      END        
