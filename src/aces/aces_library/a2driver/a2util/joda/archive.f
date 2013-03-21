      SUBROUTINE ARCHIVE( E,V1,V2,STEP,NCOMP,VEC)
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
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      Logical IsTher, IsOpn,  XYZIN, NWFINDIF
      CHARACTER*10 CRAP,CRAP2
      CHARACTER*(fnamelen) FNAME
      DIMENSION V1(NXM6),V2(NXM6*NXM6),STEP(NXM6), VEC(NOPT)
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
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



      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $        ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $        XYZTol
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
C SQUASH INTERNAL COORDINATE VECTOR IF NEEDED
C
      Write(6,*)
      Print*, "The COORD COMMON BLOCK AT ARCHIVE"
      Write(6, "(3F10.5)") (Q (I), I= 1, 3*NATOMS)
      Write(6,*)
      If (.not. xyzin) Write(6, "(3F10.5)") (R(I), I= 1, 3*NATOMS)
      Write(6,*)
C
      IF(NCOMP.EQ.1 .AND. .NOT. XYZIN) CALL SQUSH(R,NX)
C
C Open up the archive -- on the 0th cycle, we want to get rid
C of any previous file.  On subsequent cycles, we want to
C open for append.
C
      call gfname(arcfil,fname,ilength)
      If (Ncycle .eq. 0) then
       Inquire (File=fname(1:ilength), Exist=IsTher)
       If (IsTher) then
        Inquire(File=fname(1:ilength), Opened=IsOpn, Number=LuOld)
        If (IsOpn) then
         Close(LuOld, Status='DELETE')
        Else
         Open(LuArc,File=fname(1:ilength),Form='UNFORMATTED',
     &        Status='OLD')
         Close(LuArc, Status='DELETE')
        EndIf
       EndIf
       OPEN(UNIT=LUARC,FILE=fname(1:ilength),FORM='UNFORMATTED',
     &      STATUS='NEW')
      ElseIf (NCYCLE .gt. 0) THEN
       OPEN(UNIT=LUARC,FILE=fname(1:ilength),FORM='UNFORMATTED',
     &      STATUS='OLD')
       DO 10 J=1,NCYCLE
        READ(LuArc)CRAP2
        IF(CRAP2.NE.CRAP) then
         Write (LuErr,*) ' *Problem with optimization ',
     &          'archive file.'
         Call ErrEx
        EndIf
 10    Continue
      ENDIF
C
      Write(6,*) "@-ARCHIVE, the number of opt. cycles", NCYCLE
      WRITE(LuArc)CRAP,iarch,NCYCLE,NATOMS,NX,NUNIQUE,NOPT,IPRNT,INR,
     $     IVEC,
     $   IDIE,IMXSTP,ISTCRT,IVIB,ICURVY,ICONTL,IRECAL,INTTYP,IDISFD,
     $   IGRDFD,ICNTYP,ISYM,IBASIS,
     $   E,(V1(I),I=1,NXM6),(V2(I),I=1,NXM6*NXM6),
     $   (STEP(J),J=1,NXM6),
     $   (R(J),J=1,NXM6),(Q(J),J=1,NX),(NCON(I),I=1,NX),
     $   (ZSYM(I),I=1,NATOMS),
     $   (VARNAM(I),I=1,NX),(PARNAM(I),I=1,NX),
     $   (ISQUASH(J),J=1,NXM6),
     $   (IATNUM(I),I=1,NATOMS),(ATMASS(I),I=1,NATOMS),
     $   (IUNIQUE(J),J=1,NUNIQUE),
     $   (NEQ(J),J=1,NXM6),
     $   ((IEQUIV(I,J),J=1,NEQ(IUNIQUE(I))),I=1,NUNIQUE),
     $   (NOPTI(J),J=1,NOPT), (VEC(J),J=1,NOPT),FPGRP,BPGRP,PGRP,
     $   ((ORIENT(I,J),J=1,3),I=1,3)
      CLOSE(LuArc,STATUS='KEEP')
C
      RETURN
      END
