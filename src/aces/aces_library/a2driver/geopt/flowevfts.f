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
      SUBROUTINE FLOWEVFTS(SCRATCH, GRDHES, HESMOD, DIAGHES, LMBDAN, 
     &                     LMBDAP, STPMAX, MORSE, IMODE, NX, LUOUT, 
     &                     IBREAK, NOPT, QSTLST_CLIMB)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DOUBLE PRECISION LMBDAN, LMBDAP
      LOGICAL MORSE, QSTLST_CLIMB
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), HESMOD(NOPT, NOPT),
     &          DIAGHES(NOPT, NOPT)
C
      IBREAK = 0
C
      IF (.NOT.QSTLST_CLIMB) THEN
C
      DO 10 I = 1, NOPT
C
         IF (IBREAK .EQ. 1) GOTO 10
C
         IF (I .EQ. IMODE) THEN
C
C Allow here for the possibility that the user may want break symmetry
C following a non-totally symmetry mode out of a local minimum. In this
C case Lambda(P) and eigenvalues, gradient along the cooresponding 
C eigen vectors become zero. If this occurs, the code below will force
C the geometry to follow this mode.
C
            DENOM = HESMOD(I,I) - LMBDAP
C
            IF (DABS(GRDHES(I)) .LT. 1.0D-8 .AND. 
     &           DENOM .LT. 1.0D-11) THEN
               CALL ZERO(SCRATCH(NOPT + 1), NOPT)
               GRDHES(I) = 100.D0
               DENOM = 1.0D0
C
               WRITE(LUOUT, 4032) IMODE
 4032          FORMAT(T3,'@EFOL-I, Following mode ',I3,' results in ',
     &               'symmetry lowering.',/,T3,' A small step will',
     &               ' be taken ', 'along this mode.')
C        
               IBREAK = 1
               STPMAX = 0.100
               WRITE(LUOUT, 4033)NINT(STPMAX*100)
 4033          FORMAT(T3,'@EFOL-I, Step size will be ',I3,' millibohr.')
            ENDIF
C
         ELSE
C
            DO 20 J = 1, NOPT
C
               SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-GRDHES(I)
     &                           *DIAGHES(J,I)/(HESMOD(I,I) - LMBDAN)
C
 20         CONTINUE
C
         ENDIF
C
 10   CONTINUE
C
C     ENDIF (.NOT.QSTLST_CLIMB)
      ENDIF
C
      IF (MORSE) CALL DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
C
C Add in part of step which goes along the reaction coordinate.
C
      DO 30 J = 1, NOPT
         IF (QSTLST_CLIMB) DENOM = HESMOD(IMODE, IMODE) - LMBDAP
C
         SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-GRDHES(IMODE)
     &                     *DIAGHES(J,IMODE)/DENOM
C     
 30   CONTINUE
CSSS      WRITE(6,*) "The unscaled step size", (SCRATCH(J+NOPT),
CSSS     &            J=1, NOPT)
C
      RETURN
      END
