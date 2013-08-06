
      SUBROUTINE FOLOWDFLT(HESMOD, DIAGHES, SCRATCH, VEC, TS, NRORMANR,
     &                     RFA, IVEC, IMODE, NCYCLE, NX, NOPT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL TS, NRORMANR, RFA

c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end

      DIMENSION HESMOD(NOPT, NOPT), DIAGHES(NOPT, NOPT),
     &          SCRATCH(NX*NX), VEC(NOPT)

      DATA IONE /1/









































































































































































































































































































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/





      CALL GETREC(-1,'JOBARC','FIRSTCAL',IONE,ICYCLE)

      INEG = 0
      DO 1101 J=1, MAX(IVEC, 1)
         Z=100.D0
         DO 1100 I = 1, NOPT
            IF (SCRATCH(I).LT.Z) THEN
c              Assignment of eigenvector to be followed
               IF (ICYCLE.EQ.0.AND.TS) IMODE = I
               IMODE = I
               Z = SCRATCH(I)
            ENDIF
            IF (J.EQ.1.AND.HESMOD(I,I).LT.0.D0) INEG = INEG + 1
 1100    CONTINUE
         IF (ICYCLE.EQ.0.AND.TS) SCRATCH(IMODE) = 99.0D00
 1101 CONTINUE

C We follow nothing for minimas. Do this just to be on safe grounds.
      IF (NRORMANR) IMODE = 0

C On latter passes, determine overlap between Hessian eigenvectors
C and VEC (saved from the previous run). I am not sure whether this
C is the only strategy we should follow. See above for explanation.
C Ajith Perera 07/04

      IF (ICYCLE.GT.0.AND.TS) THEN
         ZOVLP = 0.d0
         DO 1833 I=1, NOPT
            SVH = xdot(NOPT,VEC,1,DIAGHES(1,I),1)
            IF (DABS(SVH).GT.ZOVLP) THEN
               ZOVLP = DABS(SVH)
               IMODE = I
            END IF
 1833    CONTINUE
         WRITE(LUOUT, 1351)IMODE, ZOVLP
 1351    FORMAT(T3,' Eigenvector ',i2,' has largest overlap with last ',
     &          'mode followed.',/,T3,' Magnitude of overlap is ',
     &           f8.5,'.')
      END IF

      IF (INEG.GT.0.AND.NRORMANR) THEN
         IDIE = IFLAGS(50)
         IF (IDIE.EQ.0) then
            Write (LUOUT, 3131)
 3131       FORMAT(' @ANLYSHES, Negative eigenvalues in Hessian.')
            Call ERREX
         ElSE IF (IDIE.EQ.2) THEN
            WRITE(LUOUT,7141)
 7141       FORMAT(T3,' Negative eigenvalues are present in the ',
     $           'Hessian.'/t3,' Mode switched to rational function ',
     $           'approximation.')
            RFA = .TRUE.
            NRORMANR = .FALSE.
            IMODE = 0
         ELSE IF (IDIE.EQ.1) THEN
            DO 2411 I=1, NOPT
               HESMOD(I,I) = DABS(HESMOD(I,I))
 2411       CONTINUE
            WRITE(LUOUT,1241)
 1241       FORMAT(T3,' Negative eigenvalues present in Hessian have ',
     &         'been converted to their absolute values.')
         END IF
      ELSE IF (INEG.GT.0.AND.RFA) THEN
            IMODE = 0
      ENDIF

c Save the eigenvector being followed for use in the next step
      DO 10 I = 1, NOPT
         VEC(I) = DIAGHES(I, IMODE)
 10   CONTINUE

C Now clean up the scratch array for future use.
      CALL PUTREC(1,'JOBARC','FIRSTCAL',IONE,NCYCLE)
      CALL ZERO(SCRATCH, NX*NX)

      RETURN
      END

