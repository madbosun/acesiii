      SUBROUTINE ROTM( IAXIS, ANG, IX, RT )
C
C     FORMS THREE DIMENSIONAL ROTATION MATRIX (IX=0 IF ANG IN RADS,
c     IX=1 IF DEG,  ADD 10 TO IX TO GET TRANSPOSE)
cend
*mdc*if cray
*      implicit none
*mdc*endif
      implicit double precision (a-h,o-z)
      dimension rt(3,3)
      integer ix
c      double precision ang, RT(3,3)
c
      integer iaxis, iaxis1, iaxis2
c      double precision sign, tang
C     
c
      CALL ZERO(RT,9)
      TANG=ANG
      IF ( mod(IX,10) .EQ. 1 ) TANG = ANG * DACOS(-1.D0)/180.D0
c
      sign = 1.0 d0
      if ( ix/10 .gt. 0 ) sign = -sign
c
c     next two axes
c     --  really mod( (iaxis-1)+1, 3 ) + 1 (clock, not mod, arithmetic)
      iaxis1 = mod(iaxis ,3)+1
      iaxis2 = mod(iaxis1,3)+1
c
      rt(iaxis,iaxis) = 1.0 d0
      RT(iaxis1, iaxis1) = COS(TANG)
      RT(iaxis2, iaxis2) = COS(TANG)
c
      rt(iaxis1,iaxis2) = sin(tang) * sign
      rt(iaxis2,iaxis1) = -rt(iaxis1,iaxis2)
c
      RETURN
      END
