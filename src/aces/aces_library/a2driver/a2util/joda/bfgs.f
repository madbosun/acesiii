
c The Broyden-Fletcher-Goldfarb-Shanno update.

      SUBROUTINE BFGS(V,H,SCRATCH,STEP,TBT)
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
     $          STEP(NXM6)

      N2=1+NXM6
      N3=1+2*NXM6
      N4=1+3*NXM6
      NX2=1+NXM6*NXM6
C
C SCRATCH(N2) = DG
C
      CALL VADD(SCRATCH(N2),V,SCRATCH(1),NXM6,-1.D0)
C
C TBT(1) =  [1/DG{dot}DX][DGxDG]
C
      CALL MODMATMUL(TBT(1),SCRATCH(N2),SCRATCH(N2),NXM6,1,NXM6,
     $               NXM6,1,NXM6)
      X0=1.D0/xdot(NXM6,SCRATCH(N2),1,STEP,1)
      CALL xscal(NXM6*NXM6,X0,TBT(1),1)
C
C SCRATCH(N3) = HDX, SCRATCH(N4) = DXH, TBT(NX2)= HDX{x}DXH
C
      CALL MODMATMUL(SCRATCH(N3),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
      CALL MODMATMUL(SCRATCH(N4),STEP,H,1,NXM6,NXM6,1,NXM6,NXM6)
      CALL MODMATMUL(TBT(NX2),SCRATCH(N3),SCRATCH(N4),NXM6,1,NXM6,
     $               NXM6,1,NXM6)
C
C Z0 = 1/DXHDX, TBT(NX2) = [1/DXHDX][HDX{x}DXH]
C
      Z0=1.D0/xdot(NXM6,STEP,1,SCRATCH(N3),1)
      CALL xscal(NXM6*NXM6,Z0,TBT(NX2),1)
C
C H(updated) = H(old) + [1/DG{dot}DX][DGxDG] - [1/DXHDX][HDX{x}DXH]
C
      CALL VADD(TBT(1),TBT(1),TBT(NX2),NXM6*NXM6,-1.D0)
      CALL VADD(H,H,TBT(1),NXM6*NXM6,1.D0)
C
C DON'T WRITE ENTIRE HESSIAN UNLESS SPECIFICALLY REQUESTED.
C FOR NOW, JUST USE AN IN-CODE PARAMETER.
C
      RETURN
      END

