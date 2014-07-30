      SUBROUTINE HUPDATE(V,H,SCRATCH,STEP,TBT,ITYPE)
C
C     ARGUMENTS PASSED ARE NEW GRADIENT, HESSIAN (TO BE READ AND
C     RETURNED UPDATED), SCRATCH ARRAY- OLD GRADIENT READ INTO FIRST
C     NXM6 POSITIONS, UPDATE TYPE
C     TBT IS JUST A SCRATCH ARRAY FOR MATRIX MULTIPLY OPERATIONS
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
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
c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)










































































































































































































































































































































































































































































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




      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION V(NXM6),H(NXM6*NXM6),SCRATCH(NX*NX),TBT(3*NXM6*NXM6),
     $   STEP(NXM6)
      LOGICAL POWEL_DONE
      CHARACTER*6 MSPORPSB
      CHARACTER*10 TRASH
      CHARACTER*(fnamelen) FNAME
      INTEGER    LuArc
      PARAMETER (LuArc  = 77)
      CHARACTER*(*) ArcFil
      PARAMETER    (ArcFil = 'OPTARC ')
      PARAMETER (TOLRNCE = 1.0D-12)
      DATA ONE /1.00D0/
C
      CALL GFNAME(ARCFIL,FNAME,ILENGTH)
      OPEN(UNIT=LUARC,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
      DO 10 I=1,NCYCLE
 10   READ(LUARC)TRASH
C
C     Most of the stuff read here is not used - hence the stupid vars.
C The previous cycle's gradients and Hessians are read into scratch and H.
C The incoming variables V, H, and STEP contain the current gradient, Hessian,
C and step size.
C
      READ(LUARC)TRASH,ICRAP,icycle,NJ,NJN,NJNK,NOPT,NJ1,NJ2,NJ3,NJ4,
     $   NJ5,NJ6,NJ7,NJ8,NJ9,NJA,NJB,NJC,NJD,NJE,NJF,NJG,E,
     &   (SCRATCH(I),I=1,NXM6),
     &   (H(I),I=1,NXM6*NXM6),(STEP(J),J=1,NXM6)
      CLOSE(UNIT=LUARC)
      IF(ITYPE.EQ.0)THEN
         WRITE(LuOut,70) icycle
 70      FORMAT(T3,' Hessian from cycle ',i2,' read and not updated.')
      ELSEIF(ITYPE.EQ.1)THEN
C
C        Powell update
C
         WRITE(LuOut,71) icycle
 71      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Powell update using ',
     &      'last two gradients and previous step.')
C
         PHI = ONE
         CALL POWEL(V, H, SCRATCH, STEP, TBT, PHI)
C
      ELSEIF(ITYPE.EQ.2)THEN
C
C        Broyden-Fletcher-Goldfarb-Shanno update
C
         WRITE(LuOut,72) icycle
 72      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' BFGS update using ',
     &      'last two gradients and previous step.')
C
         CALL BFGS(V, H, SCRATCH, STEP, TBT)

C
      ELSEIF(ITYPE.EQ.3)THEN
C
C        Murtagh-Sargent update
C
         WRITE(LuOut,73) icycle
 73      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Murtagh-Sargent update using ',
     &      'last two gradients and previous step.')
C
         POWEL_DONE = .TRUE.
         MSPORPSB   = "MSP   "
         CALL MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                POWEL_DONE, MSPORPSB)
C
      ELSEIF(ITYPE.EQ.4)THEN
C
C        (1-phi)*Powell + (phi)*Murtagh-Sargent update
C
         WRITE(LuOut,74) icycle
 74      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Bofill mixture of Powell and Murtagh-Sargent',
     &      ' update using last two gradients and previous step.')
C
         N2 = NXM6 + 1
         N3 = 2*NXM6+1
         N4 = 3*NXM6+1

         CALL MODMATMUL(SCRATCH(N2),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
         CALL DCOPY(NXM6, V, 1, SCRATCH(N3), 1)
         CALL DAXPY(NXM6, -ONE, SCRATCH(1), 1, SCRATCH(N3), 1)
         CALL VADD(SCRATCH(N4),SCRATCH(N3),SCRATCH(N2),NXM6,-1.D0)

         DDXDX  = DDOT(NXM6, STEP, 1, STEP, 1)
         DHDXDG = DDOT(NXM6, SCRATCH(N4), 1, SCRATCH(N4), 1)
         DHDXDX = DDOT(NXM6, SCRATCH(N4), 1, STEP, 1)

         IF (DABS(DHDXDX).GT.TOLRNCE.AND.
     &       DABS(DHDXDG).GT.TOLRNCE) THEN
            PHI = (DHDXDX**2)/(DDXDX*DHDXDG)
         ELSE
            PHI = ONE
         END IF

         CALL POWEL(V, H, SCRATCH, STEP, TBT, ONE-PHI)
         POWEL_DONE = .TRUE.
         MSPORPSB   = "BOFILL"
         CALL MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                POWEL_DONE, MSPORPSB)

      ELSE IF (ITYPE .EQ. 5) THEN
C
C        Powell-symmetric-Broyden update
C
         WRITE(LuOut,75) icycle
 75      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Powell-Symmetric-Broyden update using ',
     &      'last two gradients and previous step.')
C
         POWEL_DONE = .FALSE.
         MSPORPSB   = "PSB   "
         CALL MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                POWEL_DONE, MSPORPSB)

      ELSE IF (ITYPE .EQ. 6) THEN
C
C  Modified BFGS update

         WRITE(LuOut,76) icycle
 76      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' modified BFGS update using ',
     &      'last two gradients and previous step.')
C
         Print*, "The BFGS update check"
         Write(6,*) (V(I), I =1, NXM6)
         Write(6,*) (SCRATCH(I), I = 1, NXM6)
         Write(6,*) (STEP(I), I= 1, NXM6)
         CALL OUTPUT(H, 1, NXM6, 1, NXM6, NXM6,NXM6, 1 )
         Write(6,*)
CSSS         CALL MODF_BFGS(V, H, SCRATCH, STEP, TBT)
      ENDIF
C
      RETURN
      END
