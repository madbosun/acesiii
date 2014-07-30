
      SUBROUTINE POWEL(V,H,SCRATCH,STEP,TBT, PHI)
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
C
C X = 1/DX{dot}DX where DX is the step size.
C
      X=1.D0/xdot(NXM6,STEP,1,STEP,1)

      N2=NXM6+1
      N3=2*NXM6+1
      N4=3*NXM6+1
C
C SCRATCH(N2) = HDX where H is the Hessian, SCRATCH(N3) = DG (gradient
C difference), SCRATCH(N4) = (DG - HDX)
C
      CALL MODMATMUL(SCRATCH(N2),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
      CALL VADD(SCRATCH(N3),V,SCRATCH(1),NXM6,-1.D0)
      CALL VADD(SCRATCH(N4),SCRATCH(N3),SCRATCH(N2),NXM6,-1.D0)

      N2=NXM6*NXM6+1
      N3=NXM6*NXM6+N2
C
C TBT(1) = DXDX, TBT(N2) = (DG - HDX)DX
C
      CALL MODMATMUL(TBT(1),STEP,STEP,NXM6,1,NXM6,NXM6,1,NXM6)
      CALL MODMATMUL(TBT(N2),SCRATCH(N4),STEP,NXM6,1,NXM6,NXM6,1,
     &               NXM6)
C
C Take the transpose of (DG - HDX)DX, TBT(N3) = [(DG - HDX)DX]^{t}
C
      CALL MTRANSP(TBT(N2),TBT(N3),NXM6,NXM6,NXM6,NXM6)
C
C Built the (DG - HDX)dotDX and scale the DXDX, TBT(1) =
c [DXdot(DG - HDX)]DXDX
C
      dtmp = xdot(NXM6,SCRATCH(N4),1,STEP,1)
      CALL xscal(NXM6*NXM6,dtmp*X,TBT(1),1)
C
C TBT(N3) = (DG - HDX)DX + [(DG - HDX)DX]^{t} - [DXdot(DG - HDX)]DXDX/DX{dot}DX
C
      CALL VADD(TBT(N3),TBT(N3),TBT(N2),NXM6*NXM6,1.D0)
      CALL VADD(TBT(N3),TBT(N3),TBT(1),NXM6*NXM6,-1.D0)
C
C H(updated) = H(old) + [(DG - HDX)DX + [(DG - HDX)DX]^{t}]/DX{dot}DX +
C              DXdot(DG - HDX)]DXDX/[DX{dot}DX]^{2}
C
      CALL xscal(NXM6*NXM6,X,TBT(N3),1)
      CALL DAXPY(NXM6*NXM6, PHI, TBT(N3), 1, H, 1)
C
C DON'T WRITE ENTIRE HESSIAN UNLESS SPECIFICALLY REQUESTED.  FOR
C NOW, JUST USE AN IN-CODE PARAMETER.
C
      IF (IFLAGS(1).GE.10) THEN
         Print*, "-----The Powel updated Hessian-----"
         CALL HESSOUT(H,NXM6,NXM6,0)
      END IF

      RETURN
      END

