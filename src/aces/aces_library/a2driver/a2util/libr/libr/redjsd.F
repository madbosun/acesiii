      SUBROUTINE REDJSD(JSD1,JSD2,NCOORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION JSD1,JSD2
      DIMENSION JSD1(NCOORD,NCOORD),JSD2(NCOORD/2,NCOORD/2)
C
      NATOMS=NCOORD/6
      IOLD=1 
      DO 100 I=1,NATOMS
       JOLD=1
       DO 110 J=1,NATOMS
C
C x,x
C
        JSD2(3*(I-1)+1,3*(J-1)+1)=JSD1(IOLD,JOLD)+JSD1(IOLD+1,JOLD+1)+
     &            JSD1(IOLD+2,JOLD+2)
C
C x,y
C
        JSD2(3*(I-1)+1,3*(J-1)+2)=JSD1(IOLD,JOLD+1)+JSD1(IOLD+1,JOLD+3)+
     &            JSD1(IOLD+2,JOLD+4)
C
C x,z
C
        JSD2(3*(I-1)+1,3*(J-1)+3)=JSD1(IOLD,JOLD+2)
     &            +JSD1(IOLD+1,JOLD+4)+
     &            JSD1(IOLD+2,JOLD+5)
C
C y,x
C
        JSD2(3*(I-1)+2,3*(J-1)+1)=JSD1(IOLD+1,JOLD)+JSD1(IOLD+3,JOLD+1)+
     &            JSD1(IOLD+4,JOLD+2)
C
C y,y
C
        JSD2(3*(I-1)+2,3*(J-1)+2)=JSD1(IOLD+1,JOLD+1)
     &            +JSD1(IOLD+3,JOLD+3)+
     &            JSD1(IOLD+4,JOLD+4)
C
C y,z
C
        JSD2(3*(I-1)+2,3*(J-1)+3)=JSD1(IOLD+1,JOLD+2)
     &            +JSD1(IOLD+3,JOLD+4)+
     &            JSD1(IOLD+4,JOLD+5)
C
C z,x
C
        JSD2(3*(I-1)+3,3*(J-1)+1)=JSD1(IOLD+2,JOLD)
     &            +JSD1(IOLD+4,JOLD+1)+
     &            JSD1(IOLD+5,JOLD+2)
C
C z,y
C
        JSD2(3*(I-1)+3,3*(J-1)+2)=JSD1(IOLD+2,JOLD+1)
     &            +JSD1(IOLD+4,JOLD+3)+
     &            JSD1(IOLD+5,JOLD+4)
C
C z,z
C
        JSD2(3*(I-1)+3,3*(J-1)+3)=JSD1(IOLD+2,JOLD+2)
     &            +JSD1(IOLD+4,JOLD+4)+
     &            JSD1(IOLD+5,JOLD+5)
C
        JOLD=JOLD+6
C
110    CONTINUE
       IOLD=IOLD+6
100   CONTINUE
      RETURN
      END
