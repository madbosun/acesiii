      SUBROUTINE CHKECK(NRX,HBUF,ECKART)
C
C SUBROUTINE CHECKS NORMAL MODES AGAINST ECKART CONDITIONS.  PUTS
C   VALUES OF ROTATIONAL AND TRANSLATIONAL PSEUDOMOMENTA INTO A VECTOR
C   AND PASSES THIS BACK TO THE CALLING PROGRAM.  PRINTING IS AN OPTION.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end


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
 
 
      Common /Orient/ Orient(3,3)
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
      DIMENSION ECKART(2*NRX),HBUF(NRX,NRX),ZAT(3),ZA(3)
C
C CHECK EIGENVECTORS AGAINST ECKART CONDITIONS - CAN DETERMINE IF
C  IT IS A TRANSLATION, ROTATION OR VIBRATION.  USE TOLERANCES FOR
C  LINEAR AND ANGULAR MOMENTUM.  PUT VALUE INTO ECKART VECTOR
C  (TOP HALF USED FOR ROTATIONS, BOTTOM FOR TRANSLATIONS).
C
C
C FIRST, MAKE EIGENVECTOR DIMENSION = MASS * DISPLACEMENT SO THAT ROTATI
C  ECKART CONDITION CAN BE TESTED LATER.
C
 3456 DO 32 I=1,NRX
      DO 33 J=1,NRX
      HBuf(J,I)=HBuf(J,I)*DSQRT(ATMASS(1+(J-1)/3))
   33 CONTINUE
      CALL NORMAL(HBUF(1,I),NRX)
   32 CONTINUE
C
C TEST ECKART CONDITIONS.
C
      CALL ZERO(ECKART,2*NRX)
      DO 35 I=1,NRX
      CALL ZERO(ZAT,3)
C
C 1. CHECK ANGULAR MOMENTUM CRITERIA
C
       DO 36 J=1,NRX-2,3
        CALL CROSS(HBUF(J,I),Q(J),ZA,0)
        CALL VADD(ZAT,ZAT,ZA,3,1.D0)
   36  CONTINUE
       CALL ZERO(ZA,3)
       ECKART(I)=DIST(ZAT,ZA)
C
C 2. CHECK LINEAR MOMENTUM CRITERIA
C
C
        CALL ZERO(ZAT,3)
        DO 39 J=1,NRX-2,3
        CALL VADD(ZAT,ZAT,HBUF(J,I),3,1.D0)
   39  CONTINUE
       ECKART(NRX+I)=DIST(ZAT,ZA)
   38  CONTINUE
   35  CONTINUE
C
C DUMP ROTATIONAL AND VIBRATIONAL VALUES IF CALLED FOR.
C
      IF(IPRNT.GE.8)THEN
       WRITE(LUOUT,5700)
 5700  FORMAT(T3,' Linear and Angular "momenta" corresponding ',
     &'to Hessian Eigenvectors ',/,T10,'Mode',T30,'P',T50,'L')
       WRITE(LUOUT,5701)(I,ECKART(NRX+I),ECKART(I),I=1,NRX)
 5701  FORMAT(T10,I3,T23,F16.12,T43,F16.12)
      ENDIF
      RETURN
      END
