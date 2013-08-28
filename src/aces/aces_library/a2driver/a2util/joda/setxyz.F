
      SUBROUTINE SETXYZ(Q,NEWQ,ORIENT,szCOORD,IVAL,IOK,NATOMS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION NEWQ
      CHARACTER*1 szCOORD
      DIMENSION Q(*),NEWQ(*),ORIENT(3,3)
      PARAMETER (ONE=1.0D0)

      IOK=0
      IF (szCOORD.EQ.'X') THEN
         IF (Mod(IVAL,2).EQ.1) THEN
            DO J = 1,NATOMS
               Q(3*J) = NEWQ(3*J-2)
               Q(3*J-2) = NEWQ(3*J-1)
               Q(3*J-1) = NEWQ(3*J)
               ORIENT(1,2)=ONE
               ORIENT(2,3)=ONE
               ORIENT(3,1)=ONE
               IOK=1
            END DO
         END IF
      ELSE
         IF (szCOORD.EQ.'Y') THEN
            IF (Mod(IVAL/2,2).EQ.1) THEN
               DO J = 1,NATOMS
                  Q(3*J) = NEWQ(3*J-1)
                  Q(3*J-1) = NEWQ(3*J-2)
                  Q(3*J-2) = NEWQ(3*J)
                  ORIENT(3,2)=ONE
                  ORIENT(2,1)=ONE
                  ORIENT(1,3)=ONE
                  IOK=1
               END DO
            END IF
         ELSE
            IF (szCOORD.EQ.'Z') THEN
               IF (Mod(IVAL/4,2).EQ.1) THEN
                  CALL XCOPY(3*NATOMS,NEWQ,1,Q,1)
                  ORIENT(1,1)=ONE
                  ORIENT(2,2)=ONE
                  ORIENT(3,3)=ONE
                  IOK=1
               END IF
            END IF
         END IF
      END IF

      RETURN
      END

