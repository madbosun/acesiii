      SUBROUTINE ASSIGN_LABELS (IREDUNCO, TOTREDNCO, TOTNOFBND,
     &                          TOTNOFANG, TOTNOFDIH,INTLABEL)
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C Labels all the redundant coordinates according to whether they
C are bonds,angles or dihedrals
C
#include "mxatms.par"
C
      INTEGER TOTREDNCO,TOTNOFBND,TOTNOFANG,TOTNOFDIH
      DIMENSION IREDUNCO(4,TOTREDNCO)
      CHARACTER*4 INTLABEL(TOTREDNCO)
      CHARACTER*3 LABEL(MAXREDUNCO)
C
      DATA MONE/-1/, MTWO/-2/
      do i = 1, min(9,MAXREDUNCO)
         write(label(i),'(2i1)') 0, i
      end do
      do i = 10, min(99,MAXREDUNCO)
         write(label(i),'(i2)') i
      end do
      do i = 100, MAXREDUNCO
         write(label(i),'(i3)') i
      end do
C
      DO 10 I=1,TOTNOFBND
            IF (IREDUNCO(4,I).EQ.0) THEN
                 INTLABEL(I)='R'//LABEL(I)
            END IF
 10   CONTINUE
            
      DO 20 I=TOTNOFBND+1, (TOTREDNCO-TOTNOFDIH)
            IF (IREDUNCO(4,I).LE.MONE) THEN
                 INTLABEL(I)='A'//LABEL(I-TOTNOFBND)
            END IF
 20   CONTINUE
 
      DO 30 I=TOTNOFBND+TOTNOFANG+1,TOTREDNCO
            IF (IREDUNCO(4,I).GT.0) THEN
                INTLABEL(I)='D'//LABEL(I-TOTNOFBND-TOTNOFANG)
            END IF
 30   CONTINUE
  
C      We may avoid scanning 3 times the whole totrednco simply by getting
C      TOTNOFBND,TOTNOFANG,TOTNOFDIH (1-->TOTNOFBND,TOTNOFBND+1 -->...

      RETURN
      END
   
