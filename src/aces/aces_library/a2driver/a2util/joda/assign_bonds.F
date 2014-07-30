      SUBROUTINE ASSIGN_BONDS(IBNDTO, IREDUNCO, TOTREDNCO, TOTNOFBND,
     &                        NRATMS, MAXREDUNCO)
C
C This routine setup the connectivity array to define bonds. This
C array is dimensioned to (4, MAXREDUNCO).
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER TOTREDNCO, TOTNOFBND
      DIMENSION IREDUNCO(4, MAXREDUNCO), IBNDTO(NRATMS, NRATMS)
C
      DATA IZERO /0/

      TOTNOFBND = IZERO
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Print*, "The connectivity assignments in assign bonds"
      Print*, ((ibndto(j, i), j=1, nratms),I=1, nratms)
      Write(6,*)
#endif
C
      DO 10 IBNDS = 2, NRATMS
         DO 20 JBNDS = 1, IBNDS - 1
            
            IF (IBNDTO(JBNDS, IBNDS) .NE. 0) THEN
C
               TOTNOFBND = TOTNOFBND + 1
C
               IREDUNCO(1, TOTREDNCO + TOTNOFBND) = IBNDS
               IREDUNCO(2, TOTREDNCO + TOTNOFBND) = JBNDS
               IREDUNCO(3, TOTREDNCO + TOTNOFBND) = IZERO
               IREDUNCO(4, TOTREDNCO + TOTNOFBND) = IZERO
C
            ENDIF
C
 20      CONTINUE
 10   CONTINUE
C
      TOTREDNCO = TOTREDNCO + TOTNOFBND
C
      IF (TOTREDNCO .GE. MAXREDUNCO) THEN
        WRITE(6, 99) TOTREDNCO, MAXREDUNCO
 99     FORMAT(T3,'Maximum No. of redundant coordinates allowed',
     &         '/', ' exceeded','/',
     &         T3,' Required:',I3,' Current Maximum:',I3)
        CALL ERREX
      ENDIF
C
      RETURN
      END




