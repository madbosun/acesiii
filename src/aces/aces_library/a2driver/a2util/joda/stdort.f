      SUBROUTINE STDORT(VEC,SCRATCH,LEN,ICNTL)
C
C SUBROUTINE ROTATES MOLECULE TO THE PRINCIPAL AXIS ORIENTATION FOR
C  THE FULL MOLECULAR POINT GROUP.  ICNTL=0 ROTATES THE MOLECULE FROM
C  THE ACES2 ORIENTATION TO THE PAO, ANY OTHER VALUE ROTATES FROM THE
C  PAO TO ACES2 ORIENTATION.  LEN IS THE LENGTH OF THE VECTOR PASSED
C  IN.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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
      Common /Orient/ Orient(3,3)
      DIMENSION VEC(LEN),SCRATCH(LEN)
C
c
      NTIMES=LEN/3
      IOK=0
      IF(MOD(LEN,3).NE.0)THEN
       WRITE(LUOUT,133)LEN
 133   FORMAT(T3,'@STDORT-F, Vector with strange length [',i6,'] passed'
     &' to STDORT.')
       Call ErrEx
      ENDIF
      DO 30 I=1,3
       DO 31 J=1,3
       IF(ABS(ORIENT(I,J)).GT.1.D-3)IOK=IOK+1
 31    CONTINUE
 30   CONTINUE
      IF(IOK.EQ.0)THEN
       WRITE(6,5432)
 5432 FORMAT(T3,'@STDORT-W, Orientation matrix is empty.  Unit ',
     &'matrix used.')
      DO 40 I=1,3
       ORIENT(I,I)=1.D0
 40   CONTINUE
      ENDIF
      CALL ZERO (SCRATCH,LEN)
      CALL VADD (SCRATCH,SCRATCH,VEC,LEN,1.D0)
      CALL ZERO (VEC,LEN)
      IF(ICNTL.EQ.0)CALL MTRAN2(ORIENT,3)
C
      CALL MODMATMUL(VEC,ORIENT,SCRATCH,3,3,NTIMES,3,3,NTIMES)
      IF(ICNTL.EQ.0)CALL MTRAN2(ORIENT,3)
      IF(IOK.EQ.0)CALL ZERO(ORIENT,9)
      RETURN
      END
