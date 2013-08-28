      SUBROUTINE DOMORSE(SCRATCH, NOPT, NX, NATOMS, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
C     Maximum number of atoms currently allowed
      INTEGER MXATMS
      PARAMETER (MXATMS=100)
C
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     &                PARNAM(3*MXATMS)
      COMMON /COORD/ Q(3*MXATMS), R(3*MXATMS), NCON(3*MXATMS),
     &               NR(MXATMS), ISQUASH(3*MXATMS), IATNUM(MXATMS),
     &               ATMASS(MXATMS), IUNIQUE(3*MXATMS), NEQ(3*MXATMS),
     &               IEQUIV(3*MXATMS,3*MXATMS),NOPTI(3*MXATMS), NATOMSC
C
      DIMENSION SCRATCH(NX*NX)
C
      DO 1452 I = 1, NOPT
C
         IIN = 0
C     
         DO 1453 IJ = 1, NATOMS
C
C NOPTI(I) gives "SQUASHED" position - need unsquashed for connectivity
C
            IK = NOPTI(I) + 6 - 3/NOPTI(I)
            IF(IJ*3 + 1 .EQ. IK) IIN = IK
C
 1453    CONTINUE
C
         IF (IIN .NE. 0) THEN
            IAT = 1 + (IIN - 1)/3
            JAT = NCON(IIN)
            NA = IATNUM(IAT)
            NB = IATNUM(JAT)
C
            CALL MORSEA(NA, NB, A)
C     
            IF(A.EQ.0.D0)THEN
               WRITE(LUOUT,1474)ZSYM(IAT),ZSYM(JAT)
 1474          FORMAT(T3,'@EFOL-I, No Morse constant available',
     &                'for ',A5, '-',A5,'.  Default used.')
            ENDIF
C
C Scale factor calculated from internal coordinates (not
C symmetry coordinate increment)
C
            Z = SCRATCH(NOPT + I)/DSQRT(DFLOAT(NEQ(NOPTI(I)) + 1))
C
            IF (DABS(Z) .GT. 0.75D0) THEN
               WRITE(LUOUT, 7208) VARNAM(IK)
 7208          FORMAT(T3,' NR step for ',A5,' too large.  MANR',
     &                'scaling not done.')
            ELSE
               CORR = 1.D0 + 1.5D0*A*Z+2.3333333D0*(A*Z)**2
               WRITE(LUOUT, 7201)CORR,VARNAM(IK)
 7201          FORMAT(T3,' MANR scale factor for NR step is ',
     &                F5.3,' for ', A5,'.')
               SCRATCH(NOPT + I) = SCRATCH(NOPT + I)*CORR
            ENDIF
C     
         ENDIF
C
 1452 CONTINUE
C
      RETURN
      END


