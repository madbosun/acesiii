      SUBROUTINE USQUSH(V,ISIZE)
C
C "UNSQUASHES" A VECTOR OF LENGTH ISIZE.
C
      DOUBLE PRECISION V(ISIZE+6)
      DO 10 I=ISIZE,4,-1
      V(I+6)=V(I)
 10   CONTINUE
      V(8)=V(3)
      V(7)=V(2)
      V(4)=V(1)
      V(9)=0.D0
      V(6)=0.D0
      V(5)=0.D0
      DO 20 I=1,3
      V(I)=0.D0
 20   CONTINUE
      RETURN
      END
