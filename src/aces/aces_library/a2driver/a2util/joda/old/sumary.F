
      SUBROUTINE SUMARY(SCRATCH, RFAMAT, GRD, SCALE, STPMAG, VARNAM, 
     &                  ISQUASH, NOPT, NX, NXM6, MXATMS, IBREAK, 
     &                  ICONTL, LUOUT, NOPTI, LUDONE, DONFIL)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*5 VARNAM
      CHARACTER*(*) DONFIL
      LOGICAL bDonFil
C
      DIMENSION SCRATCH(NX*NX), RFAMAT(NOPT, NOPT), VARNAM(3*MXATMS), 
     &          ISQUASH(3*MXATMS), GRD(NXM6), NOPTI(3*MXATMS)
C
C Do some Printing
C
      WRITE(LUOUT, 88)STPMAG, SCALE
 88   FORMAT(T3,' Summary of Optimization Cycle: ',/,T3,
     &       ' The maximum unscaled ',
     &       'step is: ',F10.5,'.',/,T3,' Scale factor set to: ',
     &        F8.5,'.')
C
      WRITE(LUOUT, 93)
93    FORMAT(T3,' Forces are in hartree/bohr and hartree/radian.',
     &       /,T3,' Parameter values are in Angstroms and degrees.')
C
      WRITE(LUOUT, 81)
      WRITE(LuOut, 881)
      WRITE(LuOut, 81)
881   FORMAT(T3,'Parameter',T17,'dV/dR',T33,'Step',T47,'Rold',
     &       T63,'Rnew')
C
      DO 15 I=1,NOPT
C
         DO 23 J =1, I - 1
            IF (VARNAM(ISQUASH(NOPTI(I))) .EQ. 
     &          VARNAM(ISQUASH(NOPTI(J)))) GOTO 15
 23      CONTINUE

         CNF = 180.0D0/DACOS(-1.D0)
         IF (MOD(ISQUASH(NOPTI(I)),3) .EQ. 1) CNF = 0.529177249D0
         ZOP = SCRATCH(I) - SCRATCH(NOPT+I)
C
         WRITE(LUOUT, 84)VARNAM(ISQUASH(NOPTI(I))), GRD(NOPTI(I)),
     &                   SCRATCH(NOPT+I)*CNF, ZOP*CNF,
     &                   SCRATCH(I)*CNF
C
 15   CONTINUE
C
 84   FORMAT(T5,A,T10,4(F15.10,1X))
C
C Generate statistics based on gradients
C
      CALL ZERO(RFAMAT, NOPT*NOPT)
C
      DO 432 I =1, NOPT
         RFAMAT(I, 1) = GRD(NOPTI(I))
 432  CONTINUE
C
      CALL VSTAT(RFAMAT, SCRATCH(2*NOPT + 1), NOPT)
C
      WRITE(LUOUT,82) SCRATCH(2*NOPT + 2),SCRATCH(2*NOPT + 5)
 82   FORMAT(74('-'),/,T3,'Minimum force: ',F12.9,' / RMS force: ',
     &       F12.9)
 81   FORMAT(74('-'))
C
C see if the calcualtion has converged
C
      IF ( (SCRATCH(2*NOPT + 5) .LT. 10.D0**(-1*ICONTL) ) .AND.
     &     (IBREAK .EQ. 0)                                      ) THEN
C
         WRITE(6,1321)
 1321    FORMAT(80('-'))
         WRITE(6,1325)10.D0**(-1*ICONTL)
 1325    FORMAT(T3,' RMS gradient is below ',E10.5,'.')
         WRITE(6,1322)
 1322    FORMAT(T3,' Convergence criterion satisfied.  Optimization ',
     &      'completed.')
         WRITE(6,1321)

         CALL ADM

C Put file DonFil out on disk so that ACES II can terminate.
C First check to see there is any.
         INQUIRE(FILE=DonFil,EXIST=bDonFil)
         IF (bDonFil) THEN
            WRITE(*,*) '@SUMMARY: File ',DonFil,' already on disk.'
            WRITE(*,*) '          Something is very wrong.'
            INQUIRE(FILE=DonFil,OPENED=bDonFil,NUMBER=LUOLD)
            IF (bDonFil) THEN
               CLOSE(UNIT=LUOLD,STATUS='DELETE')
            ELSE
               OPEN(UNIT=LUDONE,FILE=DonFil,FORM='UNFORMATTED',
     &              STATUS='NEW')
               CLOSE(UNIT=LUDONE,STATUS='DELETE')
            END IF
         END IF

cYAU         WRITE(*,*) '@SUMMARY: Writing out completion flag to disk.'
         OPEN(UNIT=LUDONE, FILE=DonFil, FORM='UNFORMATTED',  
     &        STATUS='NEW')
         CLOSE(UNIT=LUDONE, STATUS='KEEP')

         CALL ACES_JA_FIN
         STOP

      END IF

      RETURN
      END

