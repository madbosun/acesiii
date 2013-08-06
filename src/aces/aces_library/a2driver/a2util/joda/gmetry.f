      SUBROUTINE GMETRY(We_havegeom, NOSILENT)
C
C     ROUTINE TO CALCULATE CARTESIAN COORDINATES FROM INTERNAL
C     COORDINATE REPRESENTATION.  SOME OF THIS HAS BEEN LIFTED FROM
C     PRDDO, ALTHOUGH SOME IMPROVEMENTS HAVE BEEN MADE.  UP TO 50
C     ATOMS ALLOWED.  CONNECTIVITY OF FIRST THREE MUST BE 1-2-3 IN
C     INTERNAL COORDINATE REP.
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


      LOGICAL OPTRES
      INTEGER NA(MXATMS),NB(MXATMS),NC(MXATMS)
      DIMENSION A(3,MXATMS),Q0(3*MXATMS)
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
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
C
      COMMON/RESTART2/OPTRES
C
      LOGICAL NOT_IN_BOUND, We_havegeom, NOSILENT
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
      IONE=1
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,IX)
CSSS      CALL GETREC(-1,'JOBARC', 'VIB_POPT',IONE, ipost_vib)
      If (IPrnt .ge. 20 .AND. NOSILENT) then
         Write (LuOut, *) 'GMETRY starting with R'
         Write (LuOut,'((I3,1X,f12.6))') (i,R(i),i=1,NX)
      EndIf
C
      Print*, "In Gmetry; ncycle, ix, and ipost_vib"
      Print*, ncycle, is, ipost_vib
      IF (ncycle .NE. 0. or. ix.ne.0) THEN
C
C DECOMPRESS R
C
         If (IPrnt .ge. 20 .AND. NOSILENT) Write (LuOut, '(a)')
     $      '@GMETRY-I, Decompressing R.'
         CALL USQUSH(R,NXM6)
         IF(.NOT.OPTRES .AND. NOSILENT)  WRITE(LuOut,78)
 78      FORMAT(' Updating structure...')
      ELSE
         ATOR=DACOS(-1.D0)/180.D0
         ATOB=0.529177249D0
         If (IPrnt .ge. 20 .AND. NOSILENT) Write (LuOut,*)
     $      'GMETRY: Converting to radians & bohr.',ator,atob
         DO 19 IX=8,NX-1,3
            IF(IX.NE.8)R(IX+1)=R(IX+1)*ATOR
 19         R(IX)=R(IX)*ATOR
         IF(IFLAGS(78).EQ.0)THEN
          DO 18 IX=4,NX-2,3
 18          R(IX)=R(IX)/ATOB
         ENDIF
      ENDIF
C
      If (IPrnt .ge. 20 .AND. NOSILENT) then
         Write (LuOut, *) 'GMETRY using R vector'
         Write (LuOut,'((I3,1X,f12.6))') (i,R(i),i=1,NX)
      EndIf

      CALL GEN_CART_COORD(Q0, NOSILENT)
C
C DEAL WITH CASE WHERE ANGLES OR DIHEDRALS NOT WITHIN BOUNDS HERE.
C
      CALL ZERO(Q0,3*MXATMS)
      CALL VADD(Q0,R,Q0,NX,1.D0)
      CALL XTOR(R,0)
      CALL USQUSH(R,NXM6)
      NOT_IN_BOUND = .FALSE.
      DO 107 I=1,NX
       IF(DABS(DABS(R(I))-DABS(Q0(I))).GT.1.D-4)THEN
        IF(NOSILENT) WRITE(6,8331)I,Q0(I),R(I)
 8331   FORMAT('@GMETRY-W, Internal coordinate #',i3,
     &' not within bounds:',
     &/,t3,' Value was: ',f10.5,' and has been changed to ',f10.5,'.')
        NOT_IN_BOUND = .TRUE.
       ENDIF
 107  CONTINUE
C
C Since their inception, the previous 12 lines, which were meant to
C correct poorly constructed internal coordinates, were doing more harm
C than good for geometry optimizations. I assume that this was left
C unnoticed for so many years simply because every one wrote decent
C Z-matrices. Then come ZMAT files generated automatically by MOPAC (or
C HyperChem?). Those matrices go through the 12 lines and trigger the
C "not within bound" flag, which causes joda to use its own internally-
C chosen values. The problem is that the Cartesian coordinates still
C correspond to the "bad" internal coordinates. For single point
C calculations, this is irrelevant since the internal coordinates have
C no role beyond generating the Cartesian coordinates (so why do we even
C bother for SP calcs?). However, during an optimization, we go back
C and forth between Cartesians and internals. As a result, when the
C condition that ((DABS(DABS(R(I))-DABS(Q0(I))).GT.1.D-4) is satisfied,
C the gradients are calculated for Cartesians that do NOT correspond to
C the internals that are supposed to be updated. This is a serious error
C and might partially explain why some TS searches were "meandering". To
C fix this, I created a new subroutine called GEN_CART_COORD which takes
C internal coordinates and connectivities (from ACES II common blocks)
C and generates ACES II Cartesians (nothing else). Ajith Perera, 04/2005
C
C Only one out-of-bound is all it takes to regenerate Cartesians.
      IF (NOT_IN_BOUND) CALL GEN_CART_COORD(Q0, NOSILENT)
C
      Call geomout
      RETURN
      END
