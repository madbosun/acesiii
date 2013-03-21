













































































































































































































































































































































































































































































































































      SUBROUTINE RETRIEVE(E, V1,V2,STEP,VEC)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      LOGICAL XYZIN, NWFINDIF
      CHARACTER*(fnamelen) FNAME
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
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $        ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $        XYZTol
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /INPTYP/ XYZIN,NWFINDIF
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
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      CHARACTER*10 CRAP,CRAP2
      DIMENSION V1(NXM6),V2(NXM6*NXM6),STEP(NXM6),VEC(NOPT)
      INTEGER    LuArc
      PARAMETER (LuArc  = 77)
      CHARACTER*(*) ArcFil
      PARAMETER    (ArcFil = 'OPTARC ')
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
      DATA CRAP /'***CYCLE  '/
C
      IF( ncycle .NE. 0 )WRITE(6,933)
933   FORMAT(T3,' Retrieving information from last optimization cycle.')
      CALL GFNAME(ARCFIL,FNAME,ILENGTH)
      OPEN(UNIT=LuArc,FILE=fname(1:ilength),STATUS='OLD',
     &     FORM='UNFORMATTED')
C
C NUMBER OF RECORDS IN .OPT IS ONE GREATER THAN THE NUMBER OF
C  OPTIMIZATION CYCLES.
C
      DO 101 K=1,nCYCLE
         READ(LuArc)CRAP2
         IF(CRAP2.NE.CRAP)then
            Write (LuErr,*) ' *Problem with optimization archive file.'
            Call ErrEx
         EndIf
 101  Continue
      Write(6,*)
      Print*, "Start Reading Archive file NXM6: ", NXM6
C
      READ(LuArc)CRAP2,iarchx,ICYCLE,NATOMS,NX,NUNIQUE,NOPT,IPRNT,INR,
     $     IVEC,
     $ IDIE,IMXSTP,ISTCRT,IVIB,ICURVY,ICONTL,IRECAL,INTTYP,IDISFD,
     $ IGRDFD,ICNTYP,ISYM,IBASIS,E,
     $   (V1(I),I=1,NXM6),
     $   (V2(I),I=1,NXM6*NXM6),(STEP(J),J=1,NXM6),(R(J),J=1,NXM6),
     $   (Q(J),J=1,NX),(NCON(I),I=1,NX),(ZSYM(I),I=1,NATOMS),
     &   (VARNAM(I),I=1,NX),(PARNAM(I),I=1,NX),(ISQUASH(J),J=1,NXM6),
     &   (IATNUM(I),I=1,NATOMS),(ATMASS(I),I=1,NATOMS),
     $   (IUNIQUE(J),J=1,NUNIQUE),(NEQ(J),J=1,NXM6),
     $   ((IEQUIV(I,J),J=1,NEQ(IUNIQUE(I))),I=1,NUNIQUE),
     &   (NOPTI(J),J=1,NOPT), (VEC(J),J=1,NOPT),FPGRP,BPGRP,PGRP,
     $   ((ORIENT(I,J),J=1,3),I=1,3)
         IF(CRAP2.NE.CRAP)then
            Write (LuErr,*) ' *Problem with optimization archive file.'
            Call ErrEx
         EndIf
C

C
C DECOMPRESS INTERNAL COORINATE VECTOR
C
      XYZIN = IFLAGS(68) .GT. 0
C
      IF (.NOT. XYZIN) THEN 
         DO J=NX,10,-1
            R(J)=R(J-6)
         ENDDO 
         R(8)=R(3)
         R(7)=R(2)
         R(4)=R(1)
         R(1)=0.D0
         R(2)=0.D0
         R(3)=0.D0
         R(5)=0.D0
         R(6)=0.D0
         R(9)=0.D0
       ENDIF
C

C
      CLOSE(77,STATUS='KEEP')
      RETURN
      END
