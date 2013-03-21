
      SUBROUTINE MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                   POWEL_DONE, MSPORPSB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL POWEL_DONE
      CHARACTER*6 MSPORPSB

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
     $          STEP(NXM6)
      DATA ONE /1.00D0/
C
C Murtagh-Sargent update. Note that there are conflicting formulas
C MS update. The conflicts are about the scaling factors and I
C doubt that they can cause major problems.
C
      N2 = NXM6 + 1
      N3 = 2*NXM6+1
      N4 = 3*NXM6+1

C SCRATCH(N2) = HDX where H is the Hessian
C SCRATCH(N3) = DG (gradient difference)
C SCRATCH(N4) = (DG - HDX)

      DDXDX  = DDOT(NXM6, STEP, 1, STEP, 1)
      DDXDX2 = DDXDX**2
      DHDXDG = DDOT(NXM6, SCRATCH(N4), 1, SCRATCH(N4), 1)
      DHDXDX = DDOT(NXM6, SCRATCH(N4), 1, STEP, 1)

      IF (MSPORPSB .EQ. "BOFILL") THEN
         IF (DABS(DHDXDG).GT.TOLRNCE.AND.
     &       DABS(DHDXDX).GT.TOLRNCE) THEN
            COEF    = (DHDXDX**2)/(DDXDX*DHDXDG)
            PHI_INT = ONE/DHDXDX
         ELSE
            COEF    = ONE
            PHI_INT = ONE
         END IF
      ELSE IF (MSPORPSB .EQ. "MSP   ") THEN
         PHI_INT =  ONE/DHDXDX
      ELSE IF (MSPORPSB .EQ. "PSB   ") THEN
         PHI_INT = DHDXDX/DDXDX2
      END IF

      N2 = NXM6*NXM6+1
      N3 = NXM6*NXM6+N2
C
C Form (DG - HDX)o(DG - HDX) and put it in SCRATCH(N2)
C
      CALL MODMATMUL(SCRATCH(N2), SCRATCH(N4), SCRATCH(N4), NXM6, 1,
     &               NXM6, NXM6, 1, NXM6)

      IF (MSPORPSB .EQ. "BOFILL") THEN
         CALL ZERO(TBT(N3), NXM6*NXM6)
         CALL DAXPY(NXM6*NXM6, PHI_INT, SCRATCH(N2), 1, TBT(N3), 1)
         CALL DAXPY(NXM6*NXM6, COEF, TBT(N3), 1, H, 1)
      ELSE IF (MSPORPSB .EQ. "MSP   ") THEN
         CALL DAXPY(NXM6*NXM6, PHI_INT, SCRATCH(N2), 1, H,  1)
      END IF

      IF (.NOT. POWEL_DONE .AND. (MSPORPSB .EQ. "PSB   ")) THEN
C        Form DXoDX in TBT(1)
         CALL MODMATMUL(TBT(1),STEP,STEP,NXM6,1,NXM6,NXM6,1,
     &                  NXM6)
C        Form (DG - HDX)DX in TBT(N2)
         CALL MODMATMUL(TBT(N2),SCRATCH(N4),STEP,NXM6,1,NXM6,
     &                  NXM6,1,NXM6)
C        TBT(N3) = TBT(N2)^{trnsp}
         CALL MTRANSP(TBT(N2),TBT(N3),NXM6,NXM6,NXM6,NXM6)
C        TBT(N3) = (DG - HDX)DX + ((DG - HDX)DX)^{trnsp}
         CALL VADD(TBT(N3),TBT(N3),TBT(N2),NXM6*NXM6,1.D0)
         CALL DSCAL(NXM6*NXM6, ONE/DDXDX, TBT(N3), 1)
         CALL DAXPY(NXM6*NXM6, -PHI_INT, TBT(1), 1, TBT(N3), 1)
         CALL DAXPY(NXM6*NXM6, ONE, TBT(N3), 1, H, 1)
      END IF

      IF (IFLAGS(1).GE.10) THEN
         IF (MSPORPSB .EQ. "BOFILL") THEN
             Print*, "-----The Bofill updated Hessian-----"
         ELSE IF (MSPORPSB .EQ. "PSB   ")  THEN
             Print*, "-----The PSB updated Hessian-----"
         ELSE IF (MSPORPSB .EQ. "MSP   ") THEN
             Print*, "-----The MSP updated Hessian-----"
         ENDIF
         CALL HESSOUT(H,NXM6,NXM6,0)
      END IF

      RETURN
      END

