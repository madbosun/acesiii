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
      SUBROUTINE NEWRAPH(SCRATCH, GRDHES, HESMOD, DIAGHES, MORSE, 
     &                   NOPT, NX, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "coord.com"
#include "flags.h"
#include "jodaflags.com"
C
      LOGICAL MORSE, XYZIN, NWFINDIF
C
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), DIAGHES(NOPT, NOPT),
     &          HESMOD(NOPT, NOPT)
C     
C Do the Normal Newton-Raphson update
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6, "(1x,a,I5)") "The number of degs. of freed. 
     &          at NR:", NOPT
      Write(6,*)
#endif
      DO 10 I = 1, NOPT
         DO 20 J = 1, NOPT
            IF (iFlags2(h_IFLAGS2_geom_opt) .EQ. 2 .and.
     &          HESMOD(I,I) .LT. 0.0D0) Then
                    SCRATCH(J+NOPT) = 0.0D0
            ELSE
C
                SCRATCH(J+NOPT)=SCRATCH(J+NOPT)-GRDHES(I)*
     &                          DIAGHES(J,I)/HESMOD(I,I) 
            ENDIF
 20      CONTINUE
 10   CONTINUE
C
C Notice that there is no Morse adjustments can be made for
C pure Cartesian OLPtimization (no connections), Ajith Perera,08/2008.
C
      IF (MORSE) THEN
         IF (iFlags2(h_IFLAGS2_geom_opt) .ge. 3) THEN
           CALL DOMORSEXYZ_RIC(SCRATCH, NOPT, NX, LUOUT)
         ELSE
           CALL DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
         ENDIF  
      ENDIF
C
      RETURN
      END
