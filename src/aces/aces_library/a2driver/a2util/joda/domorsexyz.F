      SUBROUTINE DOMORSEXYZ(SCRATCH, NOPT, NX, NATOMS, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      INTEGER TOTREDNCO, TOTNOFBND
      LOGICAL BONDED
C     Maximum number of atoms currently allowed
#include "mxatms.par"
C
      DIMENSION IBNDTO(MXATMS*MXATMS), IREDUNCO(4, MAXREDUNCO) 
      DIMENSION SCRATCH(NX*NX)
      LOGICAL I_UNIQUE(MAXREDUNCO)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     &                PARNAM(3*MXATMS)
      COMMON /COORD/ Q(3*MXATMS), R(3*MXATMS), NCON(3*MXATMS),
     &               NR(MXATMS), ISQUASH(3*MXATMS), IATNUM(MXATMS),
     &               ATMASS(MXATMS), IUNIQUE(3*MXATMS), NEQ(3*MXATMS),
     &               IEQUIV(3*MXATMS,3*MXATMS),NOPTI(3*MXATMS), NATOMSC
C 
      CALL GETREC(20, 'JOBARC', 'IBONDTO ', NATOMS*NATOMS, IBNDTO)
      CALL GETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
      CALL GETREC(20, 'JOBARC', 'UNIQUEDF', TOTREDNCO, I_UNIQUE)
      CALL GETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO, IREDUNCO)
C
      DO IOPT = 1,  NOPT
C
         IATM = IREDUNCO(1, IOPT)
         JATM = IREDUNCO(2, IOPT)
         KATM = IREDUNCO(3, IOPT)
C
         IF (I_UNIQUE(IOPT) .AND. KATM .EQ. 0) THEN
C
            NA = IATNUM(IATM)
            NB = IATNUM(JATM)
C
            CALL MORSEA(NA, NB, A)
C     
            IF(A.EQ.0.D0)THEN
               WRITE(LUOUT,1474)ZSYM(IATM),ZSYM(IATM)
 1474          FORMAT(T3,'@EFOL-I, No Morse constant available',
     &                'for ',A5, '-',A5,'.  Default used.')
            ENDIF
C
C Scale factor calculated from internal coordinates (not
C symmetry coordinate increment)
C
            Z = SCRATCH(NOPT + IOPT)/
     &          DSQRT(DFLOAT(NEQ(NOPTI(IOPT)) + 1))
C
            IF (DABS(Z) .GT. 0.75D0) THEN
               WRITE(LUOUT, 7208) VARNAM(IOPT)
 7208          FORMAT(T3,' NR step for ',A5,' too large.  MANR',
     &                'scaling not done.')
            ELSE
               CORR = 1.D0 + 1.5D0*A*Z+2.3333333D0*(A*Z)**2
               WRITE(LUOUT, 7201)CORR,VARNAM(IOPT)
 7201          FORMAT(T3,' MANR scale factor for NR step is ',
     &                F5.3,' for ', A5,'.')
               SCRATCH(NOPT + IOPT) = SCRATCH(NOPT + IOPT)*CORR
            ENDIF
C     
         ENDIF
C
      ENDDO 
C
      RETURN
      END


