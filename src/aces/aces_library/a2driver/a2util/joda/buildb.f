      SUBROUTINE BUILDB(B0,B,G,A,BT,SCRATCH,SCR2)
C
C CONSTRUCTS B MATRIX (R=BX) ANALYTICALLY - FORMULAS FOR DERIVATIVES
C   CALCULATED EXPLICITLY FROM VECTOR RELATIONS - GOES ON TO CALCULATE
C   G MATRIX AND A MATRIX.  THESE ARE HELD IN MEMORY. RETURNS G MATRIX
C   IN SCR2 WHEN VIBRATIONAL CALCULATION PERFORMED.  UNIT MATRIX USED
C   FOR M MATRIX UNLESS SPECIFICALLY REQUESTED TO USE REAL RECIPROCAL
C   MASS MATRIX.
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


c
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
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
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
      DIMENSION B(NXM6,NX),VECAB(3),VECBC(3),VECDC(3),A(NX,NXM6)
      DIMENSION ABXBC(3),DCXCB(3),VECCB(3),G(NXM6,NXM6),BT(NX,NXM6)
      DIMENSION B0(NX,NX),SCRATCH(NXM6*2),SCR2(NXM6,NXM6)
C
C TOLERANCE USED IN CALLS TO FILTER ROUTINE - INTENDED TO INCREASE
C   PRECISION OF CALCULATION
C
      PARAMETER (FILTOL = 1.D-14)
      IND(I,J)=I*3-2+J
c
      DET=0.D0
      Write(6,*)
      Print*, "The Cartesians, internals and connec. @-buildb" 
      Print*, "The NATOMS", NATOMS
      Write(6,"(3F10.5)"), (Q(I), I=1, NX)
      Write(6,*)
      Write(6,"(3F10.5)"), (R(I), I=1, NX)
      Write(6,"(4I5)")  (NCON(I), I=1, NX)
C
C
C FILL IN WITH DERIVATIVES OF INTERATOMIC DISTANCES WRT COORDINATE
C
      DTOR=DACOS(-1.D0)/180.D0
      DO 10 I=4,NATOMS*3-2,3
      DO 10 J=0,2
      B0(I,I+J)=(Q(I+J)-Q(IND(NCON(I),J)))/DIST(Q(I),Q(IND(NCON(I),0)))
10    B0(I,IND(NCON(I),J))=-B0(I,I+J)
C
C FILL IN WITH DERIVATIVES OF VALENCE ANGLES WRT COORDINATES
C
      DO 20 I=8,NATOMS*3-1,3
      IA=I-1
      IB=IND(NCON(I-1),0)
      IC=IND(NCON(I),0)
      RAB=DIST(Q(IA),Q(IB))
      RAC=DIST(Q(IA),Q(IC))
      RBC=DIST(Q(IB),Q(IC))
      DENOM=4.D0*RAB*RAB*RBC*RBC
      S2=RAB*RAB+RBC*RBC-RAC*RAC
      DO 20 J=0,2
      BC=(Q(IC+J)-Q(IB+J))/RBC
      BA=(Q(IA+J)-Q(IB+J))/RAB
      B0(I,IA+J)=(-4.D0*RAB*RBC*RBC*BC+S2*2.D0*RBC*BA)/
     $   (DSIN(R(I))*DENOM)
      B0(I,IC+J)=(-4.D0*RAB*RBC*BA*RAB+S2*2.D0*RAB*BC)/
     $   (DSIN(R(I))*DENOM)
20    B0(I,IB+J)=(-4.D0*RAB*RBC*(RAB*BA+RBC*BC)+S2*2.D0*(RAB*BC+RBC*BA))
     &   /(-DSIN(R(I))*DENOM)
C
C FILL IN WITH DERIVATIVES OF DIHEDRALS WRT CARTESIANS
C
      DO 30 I=12,NATOMS*3,3
      IA = I - 2
      IB=IND(NCON(I-2),0)
      IC=IND(NCON(I-1),0)
      ID=IND(NCON(I),0)
      DAB=DIST(Q(IA),Q(IB))
      DBC=DIST(Q(IB),Q(IC))
      DCD=DIST(Q(IC),Q(ID))
      CALL VEC(Q(ID),Q(IC),VECDC,1)
      CALL VEC(Q(IC),Q(IB),VECCB,1)
      CALL VEC(Q(IA),Q(IB),VECAB,1)
      CALL VEC(Q(IB),Q(IC),VECBC,1)
      S2=DSIN(DTOR*(180.D0-ANGLE(VECBC,VECAB,3)))
      S3=DSIN(DTOR*(180.D0-ANGLE(VECDC,VECCB,3)))
      C2=DCOS(DTOR*(180.D0-ANGLE(VECBC,VECAB,3)))
      C3=DCOS(DTOR*(180.D0-ANGLE(VECDC,VECCB,3)))
      CALL CROSS(VECAB,VECBC,ABXBC,0)
      CALL CROSS(VECDC,VECCB,DCXCB,0)
      DO 30 J=0,2
      B0(I,IA+J)=-ABXBC(J+1)/(S2*S2*DAB)
      B0(I,IB+J)=((DBC-DAB*C2)*ABXBC(J+1)/(DBC*DAB*S2*S2))+
     &C3*DCXCB(J+1)/(DBC*S3*S3)
      B0(I,IC+J)=((DBC-DCD*C3)*DCXCB(J+1)/(DBC*DCD*S3*S3))+
     &C2*ABXBC(J+1)/(DBC*S2*S2)
30    B0(I,ID+J)=-DCXCB(J+1)/(S3*S3*DCD)
C
C COMPRESS OUT ZERO ROWS OF B-MATRIX NOW
C
      DO 40 J=1,3*NATOMS
      B(1,J)=B0(4,J)
      IF(NATOMS.LE.2)GOTO 40
      B(2,J)=B0(7,J)
      B(3,J)=B0(8,J)
 40   CONTINUE
      DO 41 I=4,NXM6
      DO 41 J=1,NX
41    B(I,J)=B0(I+6,J)
      IF(IPRNT.GE.600)WRITE(LuOut,*)' B MATRIX '
      IF(IPRNT.GE.600)WRITE(LuOut,71)((I,J,B(I,J),J=1,NATOMS*3)
     &,I=1,NXM6)
 71   FORMAT(3('[',I3,',',I3,']',1X,F9.6,1X),'[',I3,',',I3,']',1X,
     &F9.6)
      CALL PUTREC(20,'JOBARC','BMATRIX ',IINTFP*NXM6*NX,B)
C
C SQUASH INTERNAL COORDINATE VECTOR IF NEEDED.
C
      If ( iarch .eq. 1)CALL SQUSH(R,NX)
C
C FORM G MATRIX ON THE WAY TO GETTING THE TRANSFORMATION - NOTE THAT THE
C  MASSES ARE NOT INCLUDED UNLESS THIS IS AN ENTRY FOR A NORMAL COORDINA
C  CALCULATION.
C
      CALL MTRANSP(B,BT,NXM6,NX,NXM6,NX)
C
C BUILD A MATRIX USING UNIT MATRIX FOR R. MASS MATRIX -
C  DOESN'T MATTER FOR TOTALLY SYMMETRIC INT. COORS.
C
       IF(IPRNT.GE.600)WRITE(LuOut,*)' Bt MATRIX'
       IF(IPRNT.GE.600)WRITE(LuOut,71)((I,J,
     & BT(I,J),J=1,NXM6),I=1,NX)
       CALL MODMATMUL(G,B,BT,NXM6,NX,NXM6,NXM6,NX,NXM6)
       IF(IPRNT.GE.600)WRITE(LuOut,*)' G MATRIX '
       IF(IPRNT.GE.600)WRITE(LuOut,71)((I,J,
     &G(I,J),J=1,I),I=1,NXM6)
C
C INVERT THE G MATRIX.
C
       EPS=0.0
       CALL MINV(G,NXM6,NXM6,SCRATCH,DET,EPS,0,1)
       IF(IPRNT.GE.600)  WRITE(LuOut,*)' G MATRIX DETERMINANT ',
     &DET
       IF(IPRNT.GE.600)WRITE(LuOut,*)' INVERSE G MATRIX '
       IF(IPRNT.GE.600)WRITE(LuOut,71)((I,J,
     & G(I,J),J=1,I),I=1,NXM6)
C
C FORM INVERSE TRANSFORMATION MATRIX AND EXIT (G IS NOW THE INVERSE OF G
C
C       IF(IPRNT.GE.600)WRITE(LuOut,*)' A MATRIX'
C       IF(IPRNT.GE.600)WRITE(LuOut,71)((B(I,J),J=1,NX),I=1,NXM6)
77     FORMAT((1X,F6.3))
       CALL MODMATMUL(B,G,B,NXM6,NXM6,NX,NXM6,NXM6,NX)
       CALL MTRANSP(B,A,NXM6,NX,NXM6,NX)
C
C NOTE THAT B IS OVERWRITTEN WITH A TRANSPOSE HERE - I DON'T THINK B WIL
C  BE NEEDED AGAIN
C
C Store the A matrix in JOBARC. This makes it easier to use the space
C allocated for A (Z(N3)) for other purposes. When we need the A matrix
C in Z(N3), we can restore it from the JOBARC file, Ajith Perera, 01/2005.
C
      CALL PUTREC(1, 'JOBARC', 'AMATRIX', IINTFP*NXM6*NX, A)
C
       IF(IPRNT.GE.600)WRITE(LuOut,*)' A MATRIX'
      IF(IPRNT.GE.600)WRITE(LuOut,71)((I,J,
     &A(I,J),J=1,NXM6),I=1,NX)
COBS      ENDIF
79    FORMAT((1X,F6.3))
      RETURN
      END
