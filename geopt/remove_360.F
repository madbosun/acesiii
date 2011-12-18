C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      SUBROUTINE REMOVE_360(DQ, TOTNFBND, TOTREDUNCO)
   
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER TOTNFBND, TOTREDUNCO
      DIMENSION DQ(TOTREDUNCO)
C
      TWOPI = (ATAN(DFLOAT(1))*DFLOAT(8))
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
C
#ifdef _DEBUG_LVL0
      Print*, "The Two Phi, TOTNFBND, TOTREDUNCO: ",
     &         TwOPI, TOTNFBND, TOTREDUNCO
      Write(6,*)
      Write(6,*) "Angles @Entry:"
      Write(6,*)
      Write(6,10) (DQ(I), I = (TOTNFBND + 1), TOTREDUNCO)
   10 Format (5(1X,F10.6)) 
#endif
C  
      DO i = (TOTNFBND + 1), TOTREDUNCO
C
        If(DABS(DQ(i)+TwoPi).lt.DABS(DQ(i)))
     &    DQ(i) = DQ(i) + TwoPi
        If(DABS(DQ(i)+TwoPi).lt.DABS(DQ(i)))
     &    DQ(i) = DQ(i) + TwoPi
        If(DABS(DQ(i)-TwoPi).lt.DABS(DQ(i)))
     &    DQ(i) = DQ(i) - TwoPi
        If(DABS(DQ(i)-TwoPi).lt.DABS(DQ(i)))
     &    DQ(i) = DQ(i) - TwoPi
      ENDDO
C
#ifdef _DEBUG_LVL0
      Print*, "The angles: @Exist"
      Write(6,*)
      Write(6,10) (DQ(I), I = (TOTNFBND + 1), TOTREDUNCO)
#endif
C

 
      RETURN
      END
