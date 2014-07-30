      SUBROUTINE ASSIGN_DIHLSNPL(CARTCOORD, IBNDTO, IREDUNCO, 
     &                           TOTREDNCO, TOTNOFBND, TOTNOFANG,
     &                           TOTNOFDIH, NRATMS, MAXREDUNCO,
     &                           THRESHOLD)
C
C This routine setup the connectivity array to define dihedral angles
C for nearly planar molecules. This  array is dimensioned to 
C (4, MAXREDUNCO).
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
      DIMENSION IREDUNCO(4, MAXREDUNCO), IBNDTO(NRATMS, NRATMS),
     &          CARTCOORD(3*NRATMS), TMPVEC1(3), TMPVEC2(3),
     &          TMPVEC3(3)
C
      DATA IZERO /0/
C
      NOOFNEW = IZERO
      
      DO 10 IBNDS = 1, TOTNOFBND
         IBND1 = IREDUNCO(1, IBNDS)
         IBND2 = IREDUNCO(2, IBNDS)
C
         IF (IBND2 .GT. IZERO) THEN
C
            DO 20 IANGLS = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
C
               IANG1 = IREDUNCO(1, IANGLS)
               IANG2 = IREDUNCO(2, IANGLS)
               IANG3 = IREDUNCO(3, IANGLS)
              
               IF ((IBND1 .EQ. IANG2 .AND. IBND2 .GT. IANG1 .AND. IBND2
     &             .GT. IANG3) .OR. (IBND2 .EQ. IANG2 .AND. IBND1 .GT.
     &             IANG1  .AND. IBND1 .GT. IANG3)) THEN
C     
                  ICON1 = IANG2
                  ICON2 = IANG1
                  ICON3 = IBND1
                  ICON4 = IANG3
C
                  IF (IBND1 .EQ. IANG2) ICON3 = IBND2
C
C Take care of linear angles if they are present
C     
                  IF ((IBNDTO(ICON2, ICON3)) .EQ. -99) THEN
                     ICON4 = ICON3
                     ICON3 = IANG3
                  ELSE IF ((IBNDTO(ICON3, ICON4)) .EQ. -99) THEN
                     ICON4 = ICON2
                     ICON2 = IANG1
                  ENDIF
C
C If there is an atom bonded to the ICON1 and not in the ICON2, ICON3, 
C ICON4 plane, we can define angle bends to take care of the 
C out-of-plane coordinate, no need of a dihedral angle coordinate.
C
                  CALL VEC(CARTCOORD(3*ICON3+1), CARTCOORD(3*ICON2+1), 
     &                     TMPVEC1, 0)
                  CALL VEC(CARTCOORD(3*ICON3+1), CARTCOORD(3*ICON4+1), 
     &                     TMPVEC2, 0)
                  CALL CROSS(TMPVEC1, TMPVEC2, TMPVEC3, 1)
C
                  DO 30 KBNDS = 1, TOTNOFBND
                     KBND1 = IREDUNCO(1, KBNDS) 
                     KBND2 = IREDUNCO(2, KBNDS) 
C
                     IF (KBND2 .GT. IZERO) THEN
C
                        IF (KBND1 .EQ. ICON1 .AND. KBND2 .NE. ICON2
     &                      .AND. KBND2 .NE. ICON3 .AND. KBND2 .NE.
     &                      ICON4) THEN
C
                           CALL VEC(CARTCOORD(3*ICON3 + 1), 
     &                              CARTCOORD(3*KBND2 + 1), 
     &                              TMPVEC1, 0)
                         IF (XDOT(3, TMPVEC1, 1, TMPVEC3, 1) .GT.
     &                   THRESHOLD) 
     &                   GO TO 20
                        
                        ELSE IF (KBND2 .EQ. ICON1 .AND. KBND1 .NE. 
     &                           ICON2 .AND. KBND1 .NE. ICON3 .AND.
     &                           KBND1 .NE. ICON4) THEN

                        CALL VEC(CARTCOORD(3*ICON3 + 1), 
     &                           CARTCOORD(3*KBND1 + 1), 
     &                            TMPVEC1, 0)
                        IF (XDOT(3, TMPVEC1, 1, TMPVEC3, 1) .GT. 
     &                  THRESHOLD)  
     &                  GO TO 20
C
                        ENDIF
                     ENDIF
 30               CONTINUE
C
                  NOOFNEW = NOOFNEW + 1
C     
                  IREDUNCO(1, TOTREDNCO + NOOFNEW) = ICON3
                  IREDUNCO(2, TOTREDNCO + NOOFNEW) = ICON2
                  IREDUNCO(3, TOTREDNCO + NOOFNEW) = ICON1
                  IREDUNCO(4, TOTREDNCO + NOOFNEW) = ICON4
C
               ENDIF
C            ENDIF
C
 20      CONTINUE
        ENDIF
 10   CONTINUE
C
      TOTNOFDIH = TOTNOFDIH + NOOFNEW
      TOTREDNCO = TOTREDNCO + NOOFNEW
C
      IF (TOTREDNCO .GE. MAXREDUNCO) THEN
        WRITE(6, 99) TOTREDNCO, MAXREDUNCO
 99     FORMAT(T3,'Maximum No. of redundent coordinates allowed',
     &         '/', 'exceeded','/',
     &         T3,' Required:',I3,' Current Maximum:',I3)
        CALL ERREX
      ENDIF
C
      RETURN
      END

