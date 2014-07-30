      SUBROUTINE CONVHESS(A,SCRATCH,HC,HI,AT,ID)
C
C     CONVERTS CARTESIAN HESSIAN TO INTERNAL COORDINATE HESSIAN
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




C
      CHARACTER*(fnamelen) FNAME
      LOGICAL XYZIN, NWFINDIF
      INTEGER TOTNOFBND, TOTNOFANG, TOTNOFDIH
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DOUBLE PRECISION A(NX,NXM6),HC(NX,NX),HI(NXM6,NXM6)
      DOUBLE PRECISION AT(NXM6,NX),SCRATCH(NX,NXM6)
C
      IF (ID .GT. 0) THEN
         CALL ZERO(HI, NXM6*NXM6)
         IF (XYZIN) THEN
            IF (iFlags2(5) .ge. 3) THEN
               CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
               CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
               CALL GETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
               DO IBONDS = 1, TOTNOFBND
                  HI(IBONDS, IBONDS) = 1.0D0
                END DO
                DO IANGS = TOTNOFBND + 1, TOTNOFANG + TOTNOFBND
                   HI(IANGS, IANGS) = 0.25D0
                END DO
                DO IDIHS = TOTNOFBND + TOTNOFANG + 1, TOTNOFANG 
     &                     + TOTNOFBND + TOTNOFDIH
                   HI(IDIHS, IDIHS) = 0.10D0
                END DO
                RETURN
         ELSE IF (iFlags2(5) .eq. 2) THEN
                CALL SET_2UNIT_MATRIX(HI, NX)
                RETURN
         ENDIF

        ELSE
            HI(1,1)=1.D0
            HI(2,2)=1.D0
            HI(3,3)=0.25D0
            IF(NXM6.GT.3)THEN
              DO 10 I=4,NXM6-2,3
                 HI(I,I)    =1.D0
                 HI(I+1,I+1)=0.25D0
                 HI(I+2,I+2)=0.1D0
10            CONTINUE
             ENDIF
            RETURN 
        ENDIF
      ENDIF
C
C
C The following block of code is executed only in vibrational
C frequency calculations. The purpose is to bulid a Hessian in
C internal cordinates and write it to the FCMINT file. Note the
C assumption in the transformation: the gradient is zero, hence
C no derivative of B matrix is required. This Hessian can be
C used as a starting Hessian for Transition State searches.
C In the case of frequency calculations with Cartesians let's not do
C any transformations and leave the Hessian in Cartesian
C Coordinates.  Ajith Perera 07/2003.
C
      IF (.NOT.XYZIN) THEN
         CALL MODMATMUL(SCRATCH,HC,A,NX,NX,NXM6,NX,NX,NXM6)
         CALL MTRANSP(A,AT,NX,NXM6,NX,NXM6)
         CALL MODMATMUL(HI,AT,SCRATCH,NXM6,NX,NXM6,NXM6,NX,NXM6)
      END IF
C
C Note that the Cartesian Hessian is already in the JOBARC record
C HESSINAM written by vdint for analytical Hessians and symcor for
C numerical Hessian calculations. Since there can be ordering
C issues (vmol vs ZMAT order) let's write another record. We know
C that this is in ZMAT order. 01/2006, Ajith Perera.
C
      CALL PUTREC(20,'JOBARC','CART_HES',NX*NX*IINTFP,HC)
C
      Print*, "Writing Internal or Cartesian Hessian"
      Print*, "XYZIN, NX, NXM6:", XYZIN, NX, NXM6 
      CALL OUTPUT(HC, 1, NX, 1, NX, NX, NX, 1)
      IF(ID.EQ.-1)THEN
         CALL GFNAME('FCMINT  ',FNAME,ILENGTH)
         OPEN(UNIT=50,FILE=FNAME(1:ILENGTH),FORM='FORMATTED',
     &        STATUS='UNKNOWN')
         IF (.NOT.XYZIN) THEN
            WRITE(LUOUT,112)
112         FORMAT(T3,'@CONVHESS: the FCMINT is in internals')
            WRITE(50,'((3(E18.12,1X)))')((HI(I,J),J=1,NXM6),I=1,NXM6)
C
C Write the internal Hessian to JOBARC even though
C we have it in the FCMINT file.
C
            CALL PUTREC(20,'JOBARC','INTR_HES',NXM6*NXM6*IINTFP,HI)
         ELSE
            WRITE(LUOUT,113)
113         FORMAT(T3,'@CONVHESS: Write Cartesians FCMINT file')
            WRITE(50,'((3(E18.12,1X)))')((HC(I,J),J=1,NX),I=1,NX)
         ENDIF
          CLOSE(UNIT=50,STATUS='KEEP')
      ENDIF
C
      IF(IPRNT.GE.500)THEN
       WRITE(LUOUT,110)
 110   FORMAT(T3,' @CONVHESS-I, Full Cartesian Hessian: ')
       WRITE(LuOut,90)((I,J,HC(I,J),J=1,I),I=1,NX)
       WRITE(LUOUT,111)
 111   FORMAT(T3,' @CONVHESS-I, Full internal coordinate Hessian: ')
       WRITE(LuOut,90)((I,J,HI(I,J),J=1,I),I=1,NXM6)
 90    FORMAT(3('[',I3,',',I3,']',1X,F9.6,1X),'[',I3,',',I3,']',1X,
     & F9.6)
      ENDIF
C
      RETURN
      END
