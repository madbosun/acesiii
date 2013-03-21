      SUBROUTINE XTOR(RCURVY,IPRT)
C
c TAKES CARTESIAN COORDINATES AND CONVERTS TO CURVILINEAR INTERNAL
C  COORDINATES; USED ONLY IN THE TRANSFORMATION BETWEEN RECTILINEAR
C  AND CURVILINEAR HESSIAN ELEMENTS, WHICH FOR NOW IS DONE SLOPPILY
C  (NUMERICALLY).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


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


C
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
 
 
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
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
      DIMENSION RCURVY(NX),AVECBA(3),AVECBC(3)
      DIMENSION VECAB(3),VECBC(3),VECCD(3)
      DIMENSION ABXBC(3),BCXCD(3)
      IND(I,J)=I*3-2+J
      DTOR=ACOS(-1.D0)/180.D0
C
C
C CALCULATE BOND LENGTHS
C
      DO I = 4, NATOMS*3-2, 3
         RCURVY(I) = DIST(Q(I),Q(IND(NCON(I),0)))
      END DO
C
C CALCULATE BOND ANGLES (RADIANS)
C
      DO I = 8, NATOMS*3-1, 3
         IA = I-1
         IB = IND(NCON(I-1),0)
         IC = IND(NCON(I),0)
         CALL VEC(Q(IB),Q(IA),AVECBA,1)
         CALL VEC(Q(IB),Q(IC),AVECBC,1)
         dtmp = xdot(3,AVECBA,1,AVECBC,1)
         RCURVY(I) = DACOS(dtmp)
      END DO
C
C CALCULATE DIHEDRAL ANGLES (RADIANS)
C
      DO 30 I=12,NATOMS*3,3
      IA = I - 2
      IB=IND(NCON(I-2),0)
      IC=IND(NCON(I-1),0)
      ID=IND(NCON(I),0)
      CALL VEC(Q(IC),Q(ID),VECCD,1)
      CALL VEC(Q(IA),Q(IB),VECAB,1)
      CALL VEC(Q(IB),Q(IC),VECBC,1)
      S2=DSIN(DTOR*(180.D0-ANGLE(VECBC,VECAB,3)))
      S3=DSIN(DTOR*(180.D0-ANGLE(VECCD,VECBC,3)))
      CALL CROSS(VECAB,VECBC,ABXBC,0)
      CALL CROSS(VECBC,VECCD,BCXCD,0)
C
C DETERMINE "SIGN" AND MAGNITUDE OF DIHEDRAL ANGLE NOW.
C
      AZZO = xdot(3,ABXBC,1,BCXCD,1)
      AZZO = AZZO / (S2*S3)
      IF(DABS(DABS(AZZO)-1.D0).LT.1.D-10)AZZO=SIGN(1.D0,AZZO)
C
C PUT IN FILTER TO GET RID OF DACOS(1.000000000000002) PROBLEMS
C   1/31/89:JFS
C
      IF((DABS(AZZO)-1.D0).GT.1.D-4)WRITE(6,711)AZZO
 711  FORMAT(T3,'***SOMETHING IS VERY WRONG.  VALUE OF ',F8.5,
     &' SENT TO DACOS FUNCTION IN XTOR.')
      IF(DABS(AZZO).GT.1.D0)AZZO=SIGN(1.D0,AZZO)
      A1=DACOS(AZZO)
      P = xdot(3,VECAB,1,BCXCD,1)
      P = P / (S2*S3)
      RCURVY(I)=SIGN(A1,P)
 30   CONTINUE
 
      IF(IPRT.EQ.1)THEN
       DO 90 IC2=1,NXM6
        IC=ISQUASH(IC2)
        IF(MOD(IC,3).EQ.1)THEN
         FACT=1.0D0
        ELSE
         FACT=180.0D0/DACOS(-1.0D0)
        ENDIF
        RCURVY(IC)=RCURVY(IC)*FACT
        WRITE(LUOUT,722)VARNAM(IC)(1:linblnk(VARNAM(IC))),
     &                  RCURVY(IC)
90     CONTINUE
      ENDIF
C
C SQUASH R VECTOR.
C
      CALL SQUSH(RCURVY,NX)
721   FORMAT(T3,' Internal coordinates generated by XTOR: ')
722   FORMAT(A,'=',F20.12)
      RETURN
      END
