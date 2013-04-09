      SUBROUTINE ASSIGN_ANGLS(CARTCOORD, IBNDTO, IREDUNCO, TOTREDNCO,
     &                        TOTNOFBND, TOTNOFANG, NRATMS, MAXREDUNCO, 
     &                        THRESHOLD)
C
C This routine setup the connectivity array to define bond angles. This
C array is dimensioned to (4, MAXREDUNCO).
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     DOUBLE PRECISION DFLOAT
      LOGICAL USE1, USE2, MULAOK, LOOK4_ANCHOR
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG
      DIMENSION IREDUNCO(4, MAXREDUNCO), IBNDTO(NRATMS,NRATMS),
     &          CARTCOORD(3*NRATMS)

      DATA IZERO /0/, MTWO /-2/, MONE /-1/
C     
      PI        = 180.0D0
      TOTNOFBND = TOTREDNCO
      TOTNOFANG = IZERO
      MULAOK    = .FALSE. 
      LOOK4_ANCHOR = .FALSE.
C
      DO 10 IBNDS = 2, TOTNOFBND
         IBND1 = IREDUNCO(1, IBNDS)
         IBND2 = IREDUNCO(2, IBNDS)
C         
cSSS         IF (IBND2 .GT. IZERO) THEN
cSSS            USE1 = .TRUE.
cSSS            USE2 = .TRUE.
cSSS         END IF 
C           
         DO 20 JBNDS = 1, IBNDS - 1
C
            IF (IBND2 .GT. IZERO) THEN
               USE1 = .TRUE. 
               USE2 = .TRUE.
            END IF
C
            JBND1 = IREDUNCO(1, JBNDS)
            JBND2 = IREDUNCO(2, JBNDS)
            
C I'll have to check connections from ICON1,2,3...
               
            IF (JBND2 .GT. IZERO) THEN
               IF (IBND1 .EQ. JBND1 .AND. USE1) THEN
                  ICON1 = IBND2
                  ICON2 = IBND1
                  ICON3 = JBND2
                  USE1  = MULAOK
               ELSE IF (IBND1 .EQ. JBND2 .AND. USE1) THEN
                  ICON1 = JBND1
                  ICON2 = IBND1
                  ICON3 = IBND2
                  USE1  = MULAOK
               ELSE IF (IBND2 .EQ. JBND1 .AND. USE2) THEN
                  ICON1 = IBND1
                  ICON2 = IBND2
                  ICON3 = JBND2
                  USE2  = MULAOK
               ELSE IF (IBND2 .EQ. JBND2 .AND. USE2) THEN
                  ICON1 = JBND1
                  ICON2 = IBND2
                  ICON3 = IBND1
                  USE2  = MULAOK 
               ELSE
                  GO TO 20
C
C I wish I could avoid this GoTo Statement!
C
               ENDIF
            ENDIF
C     
C For linear bonds, we need to generate two redundant coordinates for 
C for the same centers (ICON1, ICON2 and ICON3). The VLUANGLE is 
C returned in degrees from CMP_ANGLE. 
C     
            CALL CMP_ANGLE(ICON1, ICON2, ICON3, CARTCOORD, 
     &                     VLUANGLE, NRATMS)
C     
C eleminate connectivity between the edge atoms 
C
            IF (ABS(VLUANGLE - PI) .LE. THRESHOLD) THEN
                IBNDTO(ICON3, ICON1) = -99 
C
                DO 30 LBNDS = 1, NRATMS
                   IF (LBNDS .NE. ICON1 .AND. LBNDS .NE. ICON2
     &                 .AND. LBNDS .NE. ICON3) THEN

                      IF (IBNDTO(LBNDS, ICON3) .NE. IZERO) THEN
                         CALL CMP_ANGLE(ICON1, ICON3, LBNDS,
     &                                  CARTCOORD, VLUANGLE, 
     &                                  NRATMS)
                         IF (ABS(VLUANGLE - PI) .LE. THRESHOLD)
     &                   THEN
                            IBNDTO(LBNDS, ICON1) = -99
                            IBNDTO(LBNDS, ICON2) = -99
                         ENDIF
C
                      ELSE IF (IBNDTO(LBNDS, ICON1) .NE. IZERO) THEN
                         CALL CMP_ANGLE(LBNDS, ICON1, ICON3,
     &                                  CARTCOORD, VLUANGLE,
     &                                  NRATMS) 
                         IF (ABS(VLUANGLE - PI) .LE. THRESHOLD)
     &                   THEN
                            IBNDTO(ICON2, LBNDS) = -99
                            IBNDTO(ICON3, LBNDS) = -99
   
                            IF (ICON2 .GT. LBNDS)
     &                      IBNDTO(LBNDS, ICON2) = -99
                            IF (ICON3 .GT. LBNDS) 
     &                      IBNDTO(LBNDS, ICON3) = -99
                         ENDIF
C
                      ENDIF
                   ENDIF
C
 30             CONTINUE
C
C Skip assigning angles if there are connections to the middle atom.
C
                DO 40 LBNDS = 1, NRATMS
                   IF (LBNDS .NE. ICON1 .AND. LBNDS .NE. ICON2
     &                .AND. LBNDS .NE. ICON3 .AND.
     &                (IBNDTO(LBNDS, ICON2) .GT. IZERO)) GO TO 20

C Once again I would have prefered to avoid this GoTo statement!
C
 40             CONTINUE
                
                TOTNOFANG = TOTNOFANG + 1
                IREDUNCO(1, TOTNOFANG + TOTREDNCO) =
     &                   MIN0(ICON1, ICON3)
                IREDUNCO(2, TOTNOFANG + TOTREDNCO) = ICON2
                IREDUNCO(3, TOTNOFANG + TOTREDNCO) =
     &                   MAX0(ICON1, ICON3)
                IREDUNCO(4, TOTNOFANG + TOTREDNCO) = MONE
C                
C For nearly linear angles serch for another atom that is
C connected to the edge atoms. This alternative scheme was
C added with help from Peter Szalay and have the advantage
C of better handling the angles that go pass 180 during the
C optimization, Ajith Perera, 03/2011.
C
                IF (LOOK4_ANCHOR) THEN
                   DO IANCHOR = 1, NRATMS
#ifdef _DEBUG_LVL0
      Write(6, "(a,4I3)") "Check for anchor atom: ICON1", 
     &      IBNDTO(MIN(ICON1,IANCHOR), MAX(ICON1,IANCHOR)),
     &      IANCHOR, ICON1
#endif
                       IF (IBNDTO(MIN(ICON1,IANCHOR),
     &                      MAX(ICON1,IANCHOR)) .GT. 0) THEN
                            IREDUNCO(4, TOTNOFANG + TOTREDNCO)
     &                      = -IANCHOR
                            GO TO 999
                       ENDIF
#ifdef _DEBUG_LVL0
      Write(6, "(a,4I3)") "Check for anchor atom: ICON2", 
     &      IBNDTO(MIN(ICON3,IANCHOR), MAX(ICON3,IANCHOR)),
     &      IANCHOR, ICON3
#endif
                       IF (IBNDTO(MIN(ICON3,IANCHOR),
     &                     MAX(ICON3,IANCHOR)) .GT. 0) THEN
                           IREDUNCO(4, TOTNOFANG + TOTREDNCO) 
     &                     = -IANCHOR
                           GO TO 999
                       ENDIF
                   ENDDO
                ENDIF
C     
             ELSE
                ICON4 = MTWO
             ENDIF
C
 999         CONTINUE
C
             TOTNOFANG = TOTNOFANG + 1
             IREDUNCO(1, TOTNOFANG + TOTREDNCO) =
     &                MIN0(ICON1, ICON3)
             IREDUNCO(2, TOTNOFANG + TOTREDNCO) = ICON2
             IREDUNCO(3, TOTNOFANG + TOTREDNCO) = 
     &                MAX0(ICON1, ICON3)
             IREDUNCO(4, TOTNOFANG + TOTREDNCO) = MTWO
             
 20       CONTINUE
 10    CONTINUE
C
      TOTREDNCO = TOTREDNCO + TOTNOFANG
C
      IF (TOTREDNCO .GE. MAXREDUNCO) THEN
        WRITE(6, 99) TOTREDNCO, MAXREDUNCO 
 99     FORMAT(T3,'Maximum No. of redundant coordinates allowed',
     &         '/', ' exceeded','/' 
     &         T3,' Required:',I3,' Current Maximum:',I3)
        CALL ERREX
      ENDIF
C
      RETURN
      END

