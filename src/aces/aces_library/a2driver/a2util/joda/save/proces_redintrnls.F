      SUBROUTINE PROCESS_REDINTRNLS(Q, REDUNCO, IREDUNCO, IBNDTO, 
     &                              I_UNIQUE, PTGRP, TOTNOFBND, 
     &                              TOTNOFANG, TOTNOFDIH, TOTREDNCO,
     &                              NEQ, IEQUIV, NATOMS, IPLUS_MINUS,
     &                              IPRNT, EPSILON)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
C
      INTEGER TOTREDNCO, ORDER_PTGRP, TOTNOFBND, TOTNOFANG, 
     &        TOTNOFDIH
      LOGICAL SYMEQVLNT, SAME_BGN, I_UNIQUE(TOTREDNCO), SAME_END,
     &        SAME_MID, EQEV_MID, EQEV_BRDGE, SAME_TWO, ABELIAN,
     &        SAMEDIHS
      CHARACTER*4 PTGRP

      DIMENSION IREDUNCO(4, TOTREDNCO), IBNDTO(NATOMS, NATOMS),
     &          Q(3*NATOMS), NEQ(MAXREDUNCO), REDUNCO(TOTREDNCO), 
     &          IEQUIV(3*MXATMS, 3*MXATMS), 
     &          IUNIQUE_OFFSET(MAXREDUNCO), IPLUS_MINUS(3*MXATMS)
C
      CALL ZERO(IBNDTO, NATOMS*NATOMS) 
      CALL ZERO(IPLUS_MINUS, 3*MXATMS)
      ABELIAN = .FALSE. 

      IF (PTGRP .EQ. "C2v ") THEN
         ABELIAN     = .TRUE. 
         ORDER_PTGRP = 4
      ELSE IF (PTGRP .EQ. "D2h ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 8 
      ELSE IF (PTGRP .EQ. "D2  ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 4
      ELSE IF (PTGRP .EQ. "C2  ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 2
      ELSE IF (PTGRP .EQ. "C2h ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 4
      ELSE IF (PTGRP .EQ. "C s ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 2 
      ELSE IF (PTGRP .EQ. "C i ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 2 
      ELSE IF (PTGRP .EQ. "C1  ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 1  
      ENDIF 
C
C Let's get the permutation of atoms with highest order symmetry 
C operation, and the order of the highest rotation axis.
C
      IF (.NOT. ABELIAN) THEN
         CALL GETREC(20, "JOBARC", "SYMEQUIV", NATOMS, NEQ)
      ENDIF 
C 
C$$$      WRITE(6,*) "The Permu. Array"
C$$$      WRITE(6,*) (NEQ(I), I= 1, NATOMS)
C
      IF (ABELIAN) THEN
         DO IATOM = 1, NATOMS  
            DO JATOM = 1, NATOMS
               IF (IATOM .NE. JATOM) THEN
                  CALL IDNTFY_SYMEQV_ATMS(PTGRP, ORDER_PTGRP, 
     &                                    Q((IATOM - 1)*3 + 1), 
     &                                    Q((JATOM - 1)*3 + 1),
     &                                    SYMEQVLNT)
                  IF (SYMEQVLNT) THEN 
C$$$                     WRITE(6,*) IATOM, JATOM
                     IBNDTO(JATOM, IATOM) = 1
                  ENDIF 
               ENDIF
            ENDDO
         ENDDO
C
      ELSE
C
         DO IATOM = 1, NATOMS
C
            KATOM1 = NEQ(IATOM)
C
            DO JATOM = 1, IATOM - 1

               KATOM2 = NEQ(JATOM)

               IBNDTO(IATOM, KATOM1) = 1
               IBNDTO(KATOM1, IATOM) = 1
               IBNDTO(JATOM, KATOM2) = 1
               IBNDTO(KATOM2, JATOM) = 1
C
               IF (IATOM .EQ. KATOM2) THEN
                  IBNDTO(JATOM, KATOM1) = 1
                  IBNDTO(KATOM1, JATOM) = 1
               ENDIF
               IF (JATOM .EQ. KATOM1) THEN
                  IBNDTO(KATOM2, IATOM) = 1
                  IBNDTO(IATOM, KATOM2) = 1
               ENDIF
               IF (IBNDTO(JATOM, KATOM1) .EQ. 1) THEN
                  IBNDTO(KATOM2, IATOM) = 1
                  IBNDTO(IATOM, KATOM2) = 1
               ENDIF
               IF (IBNDTO(IATOM, KATOM2) .EQ. 1) THEN
                  IBNDTO(JATOM, KATOM1) = 1
                  IBNDTO(KATOM1, JATOM) = 1
               ENDIF
C
            ENDDO 
         ENDDO
C
      ENDIF
C
C Identify the symmetry unique RIC bond coordinates.
C
C$$$      WRITE(6,*) "@PROCESS REDINTERNAL"
C$$$      WRITE(6,*) IBNDTO
C$$$      WRITE(6,*) TOTNOFBND, TOTNOFANG, TOTNOFDIH
      CALL ZERO(NEQ, MAXREDUNCO)
      IUNIQUE = 0 
      DO IBONDS = 1, TOTNOFBND
         IBGN = IREDUNCO(1, IBONDS)
         IEND = IREDUNCO(2, IBONDS) 
         DO JBONDS = 1, IBONDS - 1
C$$$            IF (IPRNT .GT. 100) WRITE(6,*) IBONDS, JBONDS
            IF (IBONDS .NE. JBONDS) THEN
               SAME_BGN = .FALSE.
               SAME_END = .FALSE. 
               SAME_TWO = .FALSE.
               JBGN = IREDUNCO(1, JBONDS)
               JEND = IREDUNCO(2, JBONDS)
 
               IF (IPRNT .GT. 10) WRITE(6,*) IBGN, IEND, JBGN, 
     &                            JEND, IUNIQUE
               IF (IBGN .EQ. JBGN) SAME_BGN = .TRUE. 
               IF ((IBNDTO(JEND, IEND) .NE. 0) .AND. (SAME_BGN))
     &         THEN 
                  IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                     IOFFSET = IUNIQUE_OFFSET(JBONDS)
                  ENDIF
                  NEQ(IOFFSET) = NEQ(IOFFSET) + 1 
                  IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
C$$            WRITE(6,*) IBONDS, JBONDS, NEQ(JBONDS)
C$$            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                  GO TO 10 
               END IF
C
               IF (IEND .EQ. JEND) SAME_END = .TRUE.
               IF ((IBNDTO(JBGN, IBGN) .NE. 0) .AND. (SAME_END))
     &         THEN 
                  IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                     IOFFSET = IUNIQUE_OFFSET(JBONDS)
                  ENDIF
                  NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                  IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
C$$            WRITE(6,*) IBONDS, JBONDS, NEQ(JBONDS)
C$$            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                  GO TO 10
               ENDIF
C
               IF (IBGN .EQ. JEND) SAME_TWO = .TRUE.
               IF ((IBNDTO(JBGN, IEND) .NE. 0) .AND. (SAME_TWO))
     &         THEN
                  IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                     IOFFSET = IUNIQUE_OFFSET(JBONDS)
                  ENDIF
                  NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                  IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
C$$$            WRITE(6,*) IBONDS, JBONDS, NEQ(IUNIQUE)
C$$$            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                  GO TO 10
               ENDIF
C
               IF (IEND .EQ. JBGN) SAME_TWO = .TRUE.
               IF ((IBNDTO(JEND, IBGN) .NE. 0) .AND. (SAME_TWO))
     &         THEN
                  IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                     IOFFSET = IUNIQUE_OFFSET(JBONDS)
                  ENDIF
                  NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                  IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
C$$$            WRITE(6,*) "I am Here"
C$$$            WRITE(6,*) IBONDS, JBONDS, NEQ(IUNIQUE)
C$$$            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                  GO TO 10
               ENDIF
C
               IF ((IBNDTO(IBGN, JBGN) .NE. 0) .AND.
     &         (IBNDTO(IEND, JEND) .NE. 0) .OR.
     &         (IBNDTO(IBGN, JEND) .NE. 0) .AND.  
     &         (IBNDTO(IEND, JBGN) .NE. 0))
     &         THEN
                  IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                     IOFFSET = IUNIQUE_OFFSET(JBONDS)
                  ENDIF
                  NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                  IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
C$$            WRITE(6,*) IBONDS, JBONDS, NEQ(JBONDS)
C$$            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                  GO TO 10
               ENDIF
C 
            ENDIF
C
         END DO
C
         I_UNIQUE(IBONDS) = .TRUE. 
         IUNIQUE = IUNIQUE + 1
         IUNIQUE_OFFSET(IBONDS) = IUNIQUE
C
 10      CONTINUE
C
      END DO
C
C Identify the symmetry unique RIC angle coordinates. Also, this block 
C of code mark the bond angles that are constrained by symmetry. For
C exmple 60 degere bond angles in BF3 are symmetry constrained. 
C 
C$$$      WRITE(6,*) (IUNIQUE_OFFSET(J), J = 1, TOTNOFBND)
      DO IANGL = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
         IBGN = IREDUNCO(1, IANGL)
         IMID = IREDUNCO(2, IANGL)
         IEND = IREDUNCO(3, IANGL)
         DO JANGL = (TOTNOFBND + 1), IANGL - 1 
            SAME_MID = .FALSE.
            EQEV_MID = .FALSE. 
            JBGN = IREDUNCO(1, JANGL)
            JMID = IREDUNCO(2, JANGL)
            JEND = IREDUNCO(3, JANGL)
            IF (IANGL .NE. JANGL) THEN
C               WRITE(6,*) IANGL, JANGL, IUNIQUE_OFFSET(JANGL)
               IF (IMID .EQ. JMID)  SAME_MID = .TRUE.
                  IF (SAME_MID) THEN 
C$$$                     WRITE(6,*) "CHECKING SYM EQV ANGLES"
C$$$                     WRITE(6,*) IBGN, JBGN,IEND, JEND
C$$$                     WRITE(6,*) IBNDTO(IBGN, JBGN), IBNDTO(IEND, JEND),
C$$$     &                          IBNDTO(IEND, JBGN), IBNDTO(IBGN, JEND)
                     IF (IBGN .EQ. JBGN .AND. 
     &                   IBNDTO(IEND, JEND) .NE. 0) THEN 
C$$$                         WRITE(6,*) "I am HERE", IOFFSET
                         IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                            IOFFSET = IUNIQUE_OFFSET(JANGL)
                         ENDIF
                         NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                         IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
                         GO TO 20
                     ENDIF
                     IF (IEND .EQ. JEND .AND. 
     &                   IBNDTO(IBGN, JBGN) .NE. 0) THEN
C$$$                          WRITE(6,*) "I am HERE", IOFFSET
                         IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                            IOFFSET = IUNIQUE_OFFSET(JANGL)
                         ENDIF
                         NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                         IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
                         GO TO 20
                     ENDIF
                     IF (IBGN .EQ. JEND .AND. 
     &                   IBNDTO(IEND, JBGN) .NE. 0) THEN
C$$$                          WRITE(6,*) "I am HERE", IOFFSET
                         IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                            IOFFSET = IUNIQUE_OFFSET(JANGL)
                         ENDIF
                         NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                         IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
                         GO TO 20
                     ENDIF
                     IF (IEND .EQ. JBGN .AND. 
     &                   IBNDTO(IBGN, JEND) .NE. 0) THEN
C$$$                          WRITE(6,*) "I am HERE", IOFFSET
                         IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                            IOFFSET = IUNIQUE_OFFSET(JANGL)
                         ENDIF
                         NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                         IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
                         GO TO 20
                     ENDIF
                     IF ((IBNDTO(IBGN, JBGN) .NE. 0) .AND. 
     &                  (IBNDTO(IEND, JEND) .NE. 0)  .OR. 
     &                  (IBNDTO(IBGN, JEND) .NE. 0)  .AND.
     &                  (IBNDTO(IEND, JBGN) .NE. 0))
     &               THEN
C$$$                         WRITE(6,*) "I should be not here"
                         IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                            IOFFSET = IUNIQUE_OFFSET(JANGL)
                         ENDIF
                         NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                         IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
                         GO TO 20
                     ENDIF
                  ENDIF
C
               IF (IBNDTO(IMID, JMID) .NE. 0) EQEV_MID = .TRUE. 
               IF (EQEV_MID) THEN
                  IF ((IBNDTO(IBGN, JBGN) .NE. 0) .AND.
     &               (IBNDTO(IEND, JEND) .NE. 0)  .OR. 
     &               (IBNDTO(IBGN, JEND) .NE. 0)  .AND.
     &               (IBNDTO(IEND, JBGN) .NE. 0)) 
     &            THEN
C$$$                      WRITE(6,*) "WHY AM I HERE"
                      IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                         IOFFSET = IUNIQUE_OFFSET(JANGL)
                      ENDIF
                      NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                      IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
                      GO TO 20
                  ENDIF
C
               ENDIF
            ENDIF
C
         ENDDO
C 
         I_UNIQUE(IANGL) = .TRUE.
         IUNIQUE = IUNIQUE + 1
         IUNIQUE_OFFSET(IANGL) = IUNIQUE
C$$$         WRITE(6,*) I_UNIQUE(IANGL), IUNIQUE_OFFSET(IANGL), IANGL
C
 20      CONTINUE
C
      ENDDO 
C
C Identify the symmetry unique RIC dihedral angle coordinates. Also,
C this block of code mark the dihedrals angles that are constrained 
C by symmetry. For exmple 180 dihedral angles and dihedral angle pairs.
C
      DO IDIHS = (TOTNOFANG + TOTNOFBND + 1),  TOTREDNCO
         IBGN = IREDUNCO(1, IDIHS)
         IMD1 = IREDUNCO(2, IDIHS)
         IMD2 = IREDUNCO(3, IDIHS)
         IEND = IREDUNCO(4, IDIHS)
         DO JDIHS = (TOTNOFANG + TOTNOFBND + 1), IDIHS - 1 
C$$$            WRITE(6,*) IDIHS, JDIHS, REDUNCO(IDIHS), REDUNCO(JDIHS)
            SAMEDIHS = (DABS(REDUNCO(IDIHS) - REDUNCO(JDIHS)) 
     &                  .LE. EPSILON)
C$$$            WRITE(6,*) "THE SAME DIHEDRALS", SAMEDIHS
            IF (SAMEDIHS) THEN
                EQEV_BRDGE = .FALSE.  
                SAME_MID   = .FALSE. 
                JBGN = IREDUNCO(1, JDIHS)
                JMD1 = IREDUNCO(2, JDIHS)
                JMD2 = IREDUNCO(3, JDIHS)
                JEND = IREDUNCO(4, JDIHS)  
C
                IF (IDIHS .NE. JDIHS) THEN 
                   IF ((IMD1 .EQ. JMD1 .AND. IMD2 .EQ. JMD2) .OR. 
     &                (IMD1 .EQ. JMD2 .AND. IMD2. EQ. JMD1)) 
     &                SAME_MID  = .TRUE.
C
                   IF ((IBGN .EQ. JBGN .AND. SAME_MID .AND. 
     &                 IBNDTO(IEND, JEND) .NE. 0)) THEN
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                      GO TO 30
                   END IF 
C
                   IF ((IEND .EQ. JEND .AND. SAME_MID .AND. 
     &                 IBNDTO(IBGN, JBGN) .NE. 0)) THEN
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF
C
C If we come here, it is a plus-minus dihedral pair. Treat 
C them as if they are two independent coordianates. Note 
C that this can be improved. 
C
                       IPLUS_MINUS(IDIHS) = 1
                       GO TO 40
C
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                      GO TO 30
                   END IF 
C
                   IF ((IBNDTO(IMD1, JMD1) .NE. 0) .AND.
     &                (IBNDTO(IMD2, JMD2) .NE. 0)  .OR.
     &                (IBNDTO(IMD1, JMD2) .NE. 0)  .AND.
     &                (IBNDTO(IMD2, JMD1) .NE. 0)) 
     &                 EQEV_BRDGE = .TRUE. 
C$$$                   WRITE(6,*) "I see a EQUV BRIDGE", EQEV_BRDGE
C
                   IF (EQEV_BRDGE .AND. IEND .EQ. JEND .AND.
     &                 IBNDTO(IBGN, JBGN) .NE. 0) THEN
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                       GO TO 30
                   ENDIF
C
                   IF (EQEV_BRDGE .AND. IBGN .EQ. JBGN .AND.
     &                 IBNDTO(IEND, JEND) .NE. 0) THEN
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                       GO TO 30
                   ENDIF
   
                   IF (EQEV_BRDGE .AND. IBGN .EQ. JEND .AND.
     &                 IBNDTO(IEND, JBGN) .NE. 0) THEN
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                       GO TO 30
                   ENDIF

                   IF (EQEV_BRDGE .AND. IEND .EQ. JBGN .AND.
     &                 IBNDTO(IBGN, JEND) .NE. 0) THEN
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                       GO TO 30
                   ENDIF
C$$$                   WRITE(6,*) "I am HERE", IBNDTO(3, 6),
C$$$     &             IBNDTO(6, 3)
                   
C
                   IF ((EQEV_BRDGE .AND. IBNDTO(IBGN, JBGN)
     &                  .NE. 0) .AND.
     &                  (IBNDTO(IEND, JEND) .NE. 0) .OR. 
     &                  (IBNDTO(IBGN, JEND) .NE. 0) .AND.
     &                  (IBNDTO(IEND, JBGN) .NE. 0)) 
     &             THEN
C$$$                       WRITE(6,*) "Then I should be Here"
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
C$$$                          WRITE(6,*) "OFFSET =", IOFFSET
                       ENDIF 
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                       GO TO 30
                   END IF
C
                ENDIF
            ENDIF 
C$$            I_UNIQUE(IDIHS) = .TRUE.
C$$            IUNIQUE = IUNIQUE + 1
C
         ENDDO                
C
 40      CONTINUE
         I_UNIQUE(IDIHS) = .TRUE.
         IUNIQUE = IUNIQUE + 1
         IUNIQUE_OFFSET(IDIHS) = IUNIQUE
 30      CONTINUE
      ENDDO 
C
C$$$      WRITE(6,*) IUNIQUE_OFFSET(31), IUNIQUE_OFFSET(32),
C$$$     &           IUNIQUE_OFFSET(33), IUNIQUE_OFFSET(34)
      WRITE(6,*) I_UNIQUE
      WRITE(6,*) NEQ
C$$$      DO I = 1, IUNIQUE
C$$$         DO J =1, NEQ(I)
C$$$         WRITE(6,*) I, IEQUIV(I, J)
C$$$         ENDDO
C$$$      ENDDO
C
      RETURN
      END
