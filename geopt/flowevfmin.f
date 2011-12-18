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
      SUBROUTINE FLOWEVFMIN(SCRATCH, GRDHES, HESMOD, DIAGHES, LMBDAN,
     &                      NX, NOPT)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DOUBLE PRECISION LMBDAN
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), HESMOD(NOPT, NOPT),
     &          DIAGHES(NOPT, NOPT)
C
      DO 10 I = 1, NOPT
         DO 20 J = 1, NOPT
C
            SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-GRDHES(I)*DIAGHES(J,I)/
     &                        (HESMOD(I,I) - LMBDAN)
C
 20      CONTINUE
 10   CONTINUE 
C      
      RETURN
      END
