      SUBROUTINE NUCDIP(Q,ICHRG,NATOMS)
C
C SUBROUTINE COMPUTES NUCLEAR CONTRIBUTION TO DIPOLE MOMENT AND WRITES
C  RESULT TO DISK FILE 'NUCDIP'.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL YESNO
      DIMENSION Q(3*NATOMS),ICHRG(NATOMS),DIPNUC(3)
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
C
C COMPUTE X Y AND Z COMPONENTS.
C
      CALL ZERO(DIPNUC,3)
      DO 10 I=1,3
       DO 15 J=1,NATOMS
        DIPNUC(I)=DIPNUC(I)+FLOAT(ICHRG(J))*Q(3*(J-1)+I)
 15    CONTINUE
 10   CONTINUE
C
C CHECK TO SEE IF NUCDIP IS THERE.  IF IT IS, THEN VERY CAREFULLY
C  GET RID OF IT.
C
      INQUIRE(FILE=NDFil,EXIST=YESNO)
      IF(YESNO)THEN
       INQUIRE(FILE=NDFil,OPENED=YESNO,NUMBER=LUOLD)
       IF(YESNO)THEN
        CLOSE(LUOLD,STATUS='DELETE')
       ELSE
        OPEN(UNIT=LUNucD,FILE=NDFil,FORM='UNFORMATTED',
     &   STATUS='OLD')
        CLOSE(LUNucD,STATUS='DELETE')
       ENDIF
      ENDIF
C
C OPEN NUCDIP AND WRITE OUT THE X Y AND Z COMPONENTS.
C
      OPEN(UNIT=LUNucD,FILE=NDFil,FORM='UNFORMATTED',
     &  STATUS='NEW')
      WRITE(81)DIPNUC
      CLOSE(UNIT=LUNucD,STATUS='KEEP')
      RETURN
      END
