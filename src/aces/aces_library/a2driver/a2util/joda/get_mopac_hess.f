
      SUBROUTINE GET_MOPAC_HESS(A, SCRATCH, IMAP, HC, HI, AT, NATOMS,
     &                          NX, NXM6, NREAL, NREAL3, NREAL3M6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end








c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)
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

      DOUBLE PRECISION HC(NX*NX),SCRATCH(NX*NX), A(NX*NXM6),
     &                 AT(NXM6*NX), HI(NXM6, NXM6)
      INTEGER IMAP(NATOMS)
      COMMON /OPTCTL/ IPRNT, INR, IVEC, IDIE, ICURVY, IMXSTP, ISTCRT,
     &                IVIB, ICONTL, IRECAL, INTTYP, IDISFD, IGRDFD,
     &                ICNTYP, ISYM, IBASIS, XYZTol
      CHARACTER*(fnamelen) FNAME
      PARAMETER (IMOPAC_UNIT = 99)

      CALL GFNAME('NDDO_HES',FNAME, ILENGTH)
      OPEN(UNIT=IMOPAC_UNIT, FILE=FNAME(1:ILENGTH), FORM=
     &     "UNFORMATTED", STATUS="OLD")

      LTRIN = NREAL*3 *(NREAL*3 + 1)/2
      READ(IMOPAC_UNIT) (SCRATCH(I), I=1, LTRIN)
      CALL EXPND2(SCRATCH, HC, NREAL3)
c      Print*, "The Hessian just read from MOPAC"
CSSS      CALL OUTPUT(SCRATCH, 1, NREAL3, 1, NREAL3, NREAL3,
CSSS     &            NREAL3, 1)

C Take care of the cases that invlove dummy or ghost atoms.
      NDUMMYN = NX*NX
      IF (NREAL3.NE.NX) THEN
         NREALN  = NREAL3*NREAL3
         CALL ZERO(SCRATCH, NDUMMYN)
         CALL GETREC(1,'JOBARC','MAP2ZMAT',NATOMS,IMAP)
         IOFF = 1
         DO IATREL = 1, NATOMS
            IATMZMAT = IMAP(IATREL)
            DO IXYZ = 1, 3
               ICOL = IXYZ + (IATMZMAT-1)*3
               DO JATREL = 1, NATOMS
                  JATMZMAT = IMAP(JATREL)
                  IF (IATMZMAT .NE.0 .AND. JATMZMAT .NE. 0) THEN
                     IROW = 1 + (JATMZMAT-1)*3
                     CALL BLKCPY(HC(IOFF), 3, 1, SCRATCH, NX,
     &                           NX, IROW, ICOL)
                     IOFF = IOFF + 3
                  END IF
               END DO
            END DO
         END DO
         CALL DCOPY(NDUMMYN, SCRATCH, 1, HC, 1)
      END IF

C According to Keith Runge, the MOPAC Hessian is in mdynes/picometer.
C The conversion factor was provided by him.  Priv. comm.  12/2004
C
      CONVERT_TO_AU = (6.2414D0/27.2113957D0)*
     &                (1.0D0/0.529177249D0)**2
      CALL DSCAL(NDUMMYN, CONVERT_TO_AU, HC, 1)
C
C Convert to internal coordinates using the A matrix.
C
      CALL MODMATMUL(SCRATCH,HC,A,NX,NX,NXM6,NX,NX,NXM6)
      CALL MTRANSP(A,AT,NX,NXM6,NX,NXM6)
      CALL MODMATMUL(HI,AT,SCRATCH,NXM6,NX,NXM6,NXM6,NX,NXM6)

      IF (IPRNT.GE.10) THEN
         WRITE(LUOUT,110)
 110     FORMAT(T3,' @GET_MOPAC_HESS, Full Cartesian Hessian: ')
         CALL OUTPUT(HC, 1, NX, 1, NX, NX, NX, 1)
         Write(Luout, 111)
 111     FORMAT(T3,' @GET_MOPAC_HESS, Full internal coordinate
     & Hessian: ')
          CALL OUTPUT(HI, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
      END IF
      CLOSE(IMOPAC_UNIT)

      RETURN
      END

