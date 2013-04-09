
C GEOPT:  A Z-matrix input and geometry optimization module for
C         the ACES2 program system.
C
C By: J.F. Stanton, 7/88
C Modifications for portability by: D.E. Bernholdt 12/88
C Modifications for [name it] by: The ACES Development Team

#include "flags.h"

      PROGRAM GEOPT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "fnamelen.par"
      LOGICAL IGNORE,OPTRES,VIBRES,YESNO
      CHARACTER*(fnamelen) FNAME
C     Allocate more than enough core for the present needs
C     Estimated  by NX-6 => NX, and assume the worst!
      PARAMETER (MEMREQ = 14*9*MXATMS*MXATMS + 16*3*MXATMS +
     &           5*9*3*MXATMS + 2*3*MXATMS)
      DIMENSION Z(MEMREQ)
CJDW  Arrays for the POLYPRT routine.
CJDW  10/23/95. Should not now be needed as we are using POLYPRT0.
C      dimension xyzcor(3*mxatms),grad(3*mxatms),hess(9*mxatms*mxatms)
C      dimension iordco(  mxatms)
C
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      CHARACTER*5 ZSYM, VARNAM, PARNAM

cYAU passed down to entry for mkvmol
#include "baslen.par"
      character*(baslen) BasNam
      INTEGER TOTREDNCO

      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     1     PARNAM(3*MXATMS)
      COMMON /COORD/ Q(3*MXATMS), R(3*MXATMS), NCON(3*MXATMS),
     1     NR(MXATMS),ISQUASH(3*MXATMS),IATNUM(MXATMS),ATMASS(MXATMS),
     2     IUNIQUE(3*MXATMS),NEQ(3*MXATMS),IEQUIV(3*MXATMS,3*MXATMS),
     3     NOPTI(3*MXATMS), NATOMS
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /TOLERS/ SYMTOL,DEGTOL
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
     $     ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $     XYZTol

C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      LOGICAL  XYZIN,NWFINDIF
      Common /SymInf/ FPGrp, BPGrp, PGrp
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON/RESTART/IGNORE
      COMMON/RESTART2/OPTRES
      COMMON/RESTART3/VIBRES
      Common /Orient/ Orient(3,3)
#include "io_units.par"
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C     Use to fix orientation for C2v and D2h from one step to the next
      Integer At1, At2, XYOri
C
      integer genby(mxatms)

c ----------------------------------------------------------------------

#ifdef _TRAPFPE
c   o since JODA does not call crapsi, have it turn on FPEs
      call trapfpe
#endif /* _TRAPFPE */

      IJUNK=1
      XYZTol=1.0d-6
      SYMTOL=1.D-6
      DEGTOL=1.D-10
      ISTAT=0
      E=0.D0
      IPrnt = 0

c   o initialize the runtime environment
      call aces_init_rte

c   o evaluate and repair the health of the current file set
      call dfiles(ignore)

c   o initialize the job archive subsystem
      call aces_ja_init

c   o perform any first-run commands
      if (ignore) then
c      o JOBARC exists (somewhere)
         call putrec(1,'JOBARC','FIRSTRUN',1,0)
         call putrec(1,'JOBARC','DIRTYFLG',1,1)
      else
c      o this is the first time joda is being executed
         call putrec(1,'JOBARC','FIRSTRUN',1,1)
         call putrec(1,'JOBARC','DIRTYFLG',1,0)
         call firstrun
      end if

C IF .OPT FILE NOT THERE, READ Z-MATRIX AND GET THINGS SET UP
      CALL ENTRY(BasNam)

c      inquire(file='ALLDONE',exist=yesno)
c      if( (iflags2(138).eq.1) .and. .not. yesno)iarch=0
      call getrec(1,'JOBARC','FNDFDONE',1,iTmp)
      if ((iflags2(138).eq.1).and.(iTmp.eq.0)) iarch=0
C
C The logical argument to FETCHZ, controls whether it is the
C the very first call to the FETCHZ. In the case of user
C defined internal coordinates geometry optimizations, FETCHZ
C is called only once and this argument has no use, but in the
C case of redundent internal coordinate (RIC) optimizations, the FETCHZ
C is called repeatedly during a optimization with this argument
C set to .FALSE. During the first call it will generate RICs,
C B-matrix and A-matirx. During the next calls it will generate
C B-matrix and A-matrix that corresponds to new coordiantes using the
C RIC's generated in the first call. Ajith Perera, 08/2002
C
      IF ( iarch .NE.1) CALL FETCHZ(.TRUE., Z, MAXMEM)
C
c Marcel Nooijen 2001
c print warning if atoms are connected 3 to 1 and vibrations.
c Abort if resraman calculation
c
      if (iarch .ne. 1) then
c
         ione = 1
         ijunk = 0
         call getrec(-1, 'JOBARC', '12SWITCH', ione, ijunk)
         if (ijunk .eq. 1 .and. (iflags2(3) .eq. 1
     $        .or. ivib .ge. 1)) then
c
            write(6,*) ' *** Warning *** '
            write(6,*)
            WRITE(6,*)' Please reconstruct ZMAT file '
            WRITE(6,*)' Third atom should connect to second atom '
            WRITE(6,*)' Interchange atoms 1 and 2 in ZMAT '
            WRITE(6,*)' Replace all 1->2 and 2->1 '
            WRITE(6,*)' in third through last lines of ZMAT '
            WRITE(6,*)' Our apologies for the inconvenience'
            write(6,*)
            write(6,*)' Necessary for RESRAMAN calculations!!'
            write(6,*)' Erroneous representation of'
            write(6,*)' Normal Coordinates in internal representation'
            write(6,*)
            if (iflags2(3) .eq. 1) call errex
         endif
      endif
C
      IONE=1
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,IFNDIF)
      CALL GETREC( 0,'JOBARC','PASS1   ',IONE,IX)
      NWFINDIF=(IONE.GT.0)
      IF (NWFINDIF.AND.(.NOT.IGNORE)) THEN
        VIBRES=.true.
        IONE=1
        CALL PUTREC(20,'JOBARC','VIBRES  ',IONE,IONE)
      ELSE
        VIBRES=.false.
      ENDIF
      IONE=1
C
C Allocate memory, Note that in the case of Cartesian Optimizations
C the number of internal coordinates is not always correspond to
C NX-6, where NX=3*NATOMS. The fetchz.F dependency, gen_rednintrnls.F
C create the "REDNCORD" record. Ajith Perera, 09/2002.
C
      CALL GETREC(-1, 'JOBARC','REDNCORD', 1, TOTREDNCO)
      NXM6 = TOTREDNCO
      IF (iFlags(68).eq.0 .OR. NWFINDIF) NXM6 = MAX(NX-6, 1)
C
C FOR DUMMY B-MATRIX (AND CARTESIAN HESSIAN AFTER BUIL
    2 N1 = (NX * NX) + 1
C FOR B MATRIX (OCCUPIED BY A(TRSPSE) AFTER BUILDB.
      N2 = (NX * NX) + N1
C FOR G MATRIX (AND INTERNAL COORD HESSIAN AFT
      N3 = (NX * NX) + N2
C FOR A MATRIX
      N4 = (NX * NXM6) + N3
C FOR B TRANSPOSE
      N5 = (NX * NXM6) + N4
C FOR SCRATCH ARRAY IN BUILDB
      N6 = (NXM6 * 2) + N5
C CARTESIAN GRADIENT VECTOR
      N7 = NX + N6
C INTERNAL COORDINATE GRADIENT VECTOR
      N8 = NXM6 + N7
C INTERNAL COORDINATE HESSIAN FOR ARCHIVING
      N9 = (NXM6 * NXM6) + N8
C FOR STEP IN UPDATE
      N10 = NXM6 + N9
C FOR SCRATCH ARRAYS IN UPDATE
      N11 = ((NXM6 * NXM6) * 3) + N10
C FOR SCRATCH ARRAYS IN SYMMETRY
      N12 = (3 * NX) + N11
C    "    "       "     "    "
      N13 = NX + N12
C    "    "       "     "    "
      N14 = NX + N13
C FOR CONSTRAINED HESSIAN IN NEWRAP
      N15 = (NOPT * NOPT) + N14
C     FOR CONSTRAINED GRADIENT VECTOR IN NEWRAP
      N16 = NOPT + N15
C     CONSTRAINED ORTHOG. BASIS GRADIENT & HESSIAN IN EFOL
      N17 = NOPT + N16
      N18 = NOPT*NOPT + N17
C     Eigenvector followed on previous iter. (EFOL)
      N19 = NOPT + N18
C     Scratch matrices for EFOL - Used for augmented Hessian in RFA
C      and TS searches.  Dimension must be NOPT+1 for RFA.
      N20 = (NOPT+1)*(NOPT+1) + N19
      N21 = (NOPT+1)*(NOPT+1) + N20
      N22 = NX + N21
      N23 = NX + N22
      N24 = NX + N23
      N25 = NX*NX + N24
C
C A Bug fix (Only occured when extended to do Raman intensiies)
C Ajith 08/98, Note the generous allocation.
C
      N26 = NX*NX + N25
C
C Marcel Need a bigger scracth space.
cSSS N27 = NXM6*NXM6 + N26
C
      N27 = NX*NX + N26
C
C Allocate room for polarizabiliity derivatives Ajith & John, 08/98
C
      N28 = 9*NX + N27
      N29 = 9*NX + N28
      N30 =   NX + N29
      N31 =   NX + N30

      MEMTOP = N31
      IF(MEMTOP.GT.MEMREQ)THEN
        WRITE(LUOUT,5010)IINTFP*MEMTOP,IINTFP*MEMREQ
 5010   FORMAT(T3,'@GEOPT-F, Not enough core.',/,
     &       T3,'Words required:',I9,' Words available:',I9)
        CALL ERREX
      ENDIF

C     Dump memory information.
      IF (IPRNT .GE. 8) THEN
        WRITE (LUOUT, 9100)  MEMTOP,MEMREQ
 9100   FORMAT (T3,'@GEOPT-I, Words available: ',I9,
     &       ' Words Used ',I9,'.')
        IF(IPRNT.GT.900)THEN
          WRITE(LUOUT,9101)N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,
     &         N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,
     &         N23,N24,N25,N26,N27,N28,N29,N30,N31
 9101     FORMAT(T3,'@GEOPT-I, Offsets into core vector: ',
     &         (/,9I6))
        ENDIF
      ENDIF
C
C Gets information from OPTARC file for subsequent steps (after
C the first call)
C
      IF ( iarch .EQ.1) THEN
        CALL RETRIEVE( E, Z(N4), Z(N8), Z(1), Z(N18))
        CALL ZERO(Z, N11)
      ELSE
C
C CALL SYMMETRY PACKAGE AND NEW FINDIF SYMMETRY ROUTINES.
C
        IF (.NOT. XYZIN) CALL GMETRY
        IF (.NOT. XYZIN .OR. NWFINDIF)  CALL SYMMETRY(Z(N11), Z(N12),
     &       Z(N13))


CJDW 4/1/97. Write COMPSYOP with maximum possible length at beginning
C            of finite difference calculation.
C            Also COMPPERM, COMPCLSS.
        IF(.NOT.NWFINDIF)THEN
         CALL ZERO(Z,8*9)
         CALL PUTREC(20,'JOBARC','COMPSYOP',8*9*IINTFP,Z)
         CALL ZERO(Z,8*NATOMS)
         CALL PUTREC(20,'JOBARC','COMPPERM',8*NATOMS,Z)
         CALL ZERO(Z,8)
         CALL PUTREC(20,'JOBARC','COMPCLSS',8,Z)
        ENDIF
C
        IOFF=1+120*9+NATOMS*9
        MLTINT=IINTFP*(MEMREQ-120)
        CALL SYMDRV(Z(IOFF),Z,MLTINT,IOFF,FPGRP,'FULL')
        ITMP=IPRNT
        IF(PGRP.EQ.FPGRP)IPRNT=0
        CALL SYMDRV(Z(IOFF),Z,MLTINT,IOFF,PGRP,'COMP')
        IPRNT=ITMP
        CALL ZERO(Z,MEMREQ)
C
c SB : The if(.not.nwfindif) block used to be inside the if(natoms.gt.1)
c block.  I moved it outside so that the call to geomout is made
c even in the case of an atom.
C
        If (NAtoms .gt. 1) THEN
C
C        Call the Z-matrix analyzer - requires A matrix from BUILDB.
C
          IF (.NOT. XYZIN .OR. NWFINDIF) THEN
             CALL BUILDB(Z(1),Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),Z(N10))
             IF(.NOT.NWFINDIF)Call Analyze (Z(N10), Z(N3))
          ELSE
C
C Note, this is part is executed only in the first call to joda, and
C no need to re-run since it had alredy been run with the argument set
C to .TRUE, and the following records for the starting  geometry is in
C the JOBARC. Ajith Perera 08/2002.
C
             CALL GETREC(20, 'JOBARC','REDNCORD', 1, TOTREDNCO)
             LENGTH_GMAT=3*NATOMS*TOTREDNCO
             CALL GETREC(20, 'JOBARC', 'GMATRIX ',LENGTH_GMAT*IINTFP,
     &                   Z(N3))
CSSS             WRITE(6,*) IFLAGS2(
CSSS     &           h_IFLAGS2_geom_opt), NWFINDIF
             IF (.NOT. NWFINDIF .AND. .NOT. IFLAGS2(
     &           h_IFLAGS2_geom_opt).EQ.0) CALL ANALYZE(Z(N10), Z(N3))
          ENDIF
        endif
C
C Print out the coordinates and distance matrix
C
        IF(.NOT.NWFINDIF)THEN
          Call GeomOut
          if (natoms.gt.1) Call ADM
        EndIf
C
C IF THIS IS A FINITE DIFFERENCE CALCULATION AND ALL POINTS HAVE
C BEEN RUN, THEN CALL VIBRATIONAL ANALYSIS CODE
C
C Make the appropriate input deck
C
        IF (INTTYP .EQ. 0) THEN
cYAU          Call MkAbIn (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
cYAU     &         GenBy, Z(N11), IStat)
           write(*,*) '@GEOPT: MkAbin was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
CJDW 9/20/96. Allow for Seward (currently must make MOL file).
        ElseIf (IntTyp .eq. 1 .or. IntTyp .eq. 2
     &                        .or. IntTyp .eq. 4) then
          Call MkVMOL (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &         GenBy, Z(N11), IStat, BasNam)
        ElseIf (IntTyp .eq. 3) then
cYAU          Call MkCadp (PGrp, NAtoms, NUnique, ZSym, IAtNum, Q,
cYAU     &         GenBy, IStat)
           write(*,*) '@GEOPT: MkCadp was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
        ENDIF
        If (Mod(IStat,2) .eq. 1) then
           Call ErrEx
        ElseIf (IStat .eq. 8 .OR. IStat .eq. 10) then
          IStat = 0
c            Write (LuErr, 9810)
C           Eventually we will want to do more things here...
C           Probably want to change to program to chain to
C           Probably want to write out a CON file or similar
        EndIF
C
C NOW WRITE COORDINATE INFORMATION TO INPUT FOR INTS.  THIS IS CRUDE,
C BUT THINGS HAVE TO BE DONE THIS WAY FOR NOW. (OBSOLETE LONG AGO).
C
        CALL ZERO(Z, N13)
C        CALL PUTGEOM(3)
        itmp=iflags(18)
        iflags(18)=mod(itmp,100)
        IF (NOpt .ne. 0 .AND. NAtoms .gt. 1 .or. IVib .eq. 1 .or.
     &       IVIB.EQ.2
     &       .OR. IFLAGS(18).EQ.3
     &       .OR. IFLAGS(18).EQ.4
     &       .OR. IFLAGS(18).EQ.5
     &       .OR. IFLAGS(18).EQ.8
     &       .OR. IFLAGS(18).EQ.9
     &       .OR. IFLAGS(18).EQ.10
     &       .OR. IFLAGS2(3).NE.0
     &       ) then
          CALL ARCHIVE(E, Z(N7), Z(N8), Z(NXM6 + 1), 1,
     &                 Z(N18))
        ElseIf (NOpt .eq. 0) then
          Write (LuOut, 9820) 'single-point'
        ElseIf (NAtoms .eq. 1) then
          Write (LuOut, 9820) 'atomic'
        Endif
        iflags(18)=itmp
c
c         inquire(file='ALLDONE',exist=yesno)
c         if(.not.yesno)then
        if(.not. nwfindif .or. nwfindif .and. ifndif .eq. 1)then
          IF(IPRNT.EQ.999)WRITE(LUOUT,8888)(I,Z(I),I=1,N27)
          CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)
          CALL ACES_JA_FIN
          Stop
        endif
      END IF
 9820 Format (' @GEOPT-W, Archive file not created for ',A,
     $     ' calculation.')
C
C CALL ROUTINE TO BUILD B MATRIX.
C
      IF(.NOT.XYZIN.AND.IFLAGS(68).NE.1)THEN
        CALL BUILDB(Z(1), Z(N1), Z(N2), Z(N3), Z(N4), Z(N5), Z(N10))
      ELSE
C
C This block is executed for each iteration during a optimiation and
C requires a call to FETCHZ with argument set to .FALSE. to create
C B and A matrices for the new coordinates using the IRC's genereated
C in the first call to FETCHZ. Ajith Perera 08/2002
C
        CALL FETCHZ(.FALSE., Z, MAXMEM)
        CALL GETREC(20, 'JOBARC','REDNCORD', 1, TOTREDNCO)
        LENGTH_GMAT = 3*NATOMS*TOTREDNCO
        CALL GETREC(20, 'JOBARC', 'GMATRIX ',LENGTH_GMAT*IINTFP,
     &              Z(N3))
        CALL GETREC(20, 'JOBARC', 'BMATRIXT', LENGTH_GMAT*IINTFP,
     &              Z(N4))
      ENDIF
C
C IF THIS IS A VIBRATIONAL CALCULATION - DO ALL THE STUFF HERE AND
C THEN STOP.
C
      itmp=iflags(18)
      iflags(18)=mod(itmp,100)
      IF(IVIB.GE.1.OR.IFLAGS2(3).NE.0)THEN
C
CJDW  1/ 5/95 Write polyrate data for Rozeanne Steckler. This almost
C             certainly will not work if there are dummy atoms present.
C             In addition, it should only be called when we do analytical
C             hessians, as in the finite difference case it will not get
C             the coordinates, energy, and gradient for the reference
C             geometry. Also, POLYRATE must be ON.
CJDW 10/23/95 Change from POLYPRT to POLYPRT0 in standard code.
C
        IF(IFLAGS(54).EQ.1.AND.IFLAGS2(113).EQ.1)THEN
          CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                  .TRUE.,.TRUE.,.TRUE.,.TRUE.)
        ENDIF
C
        IHES=0
        IF(IVIB.EQ.3.OR.IFLAGS2(3).NE.0)IHES=1
        Call READGH (0,IHES,NATOMS,Z(N6),Z(1),Z(N1),Z(N2),IAVGRD,
     &       IAVHES)
        IF(ISTAT.NE.0)THEN
          WRITE(LUOUT,3432)
 3432     FORMAT(T3,'@GEOPT-F, Problem with FCM or GRD file.')
          CALL ERREX
        ENDIF
C
C PROJECT FORCE CONSTANT MATRIX AND DIPOLE DERIVATIVE MATRIX ONTO
C TOTALLY SYMMETRIC SUBSPACE FOR VIB=FINDIF CALCULATIONS
C
        IF(IVIB.EQ.3)THEN
          CALL GETREC(20,'JOBARC','FULLORDR',IONE,IORDER)
          CALL GETREC(20,'JOBARC','FULLNIRR',IONE,NIRREP)
          I000=N4
          I00A=I000+9*NATOMS*NATOMS
          I010=I00A+9*NATOMS*NATOMS
          I020=I010+9*NATOMS*NATOMS
          I030=I020+IORDER*9
          I040=I030+NATOMS*IORDER
          I050=I040+NATOMS
          I060=I050+NATOMS
          CALL SCOPY(9*NATOMS*NATOMS,Z(1),1,Z(I00A),1)
          CALL PROJFCM(NATOMS,NIRREP,IORDER,Z(I00A),
     &         Z(I010),Z(I020),Z(I030),Z(I040),Z(I050),
     &         Z(I060))
          CALL SCOPY(9*NATOMS*NATOMS,Z(I00A),1,Z(1),1)
          CALL PROJDPD(NATOMS,NIRREP,IORDER,ATMASS,Z(I00A),
     &         Z(I010),Z(I020),Z(I030),Z(I040),Z(I050),
     &         Z(I060))
          I000=N4
          I010=I000+27*NATOMS
          I020=I010+27*NATOMS
          I030=I020+IORDER*9
          I040=I030+NATOMS*IORDER
          I050=I040+NATOMS
          I060=I050+NATOMS
          I070=I060+84*NATOMS
          IF(I070.GT.MEMREQ)THEN
             WRITE(LUOUT,5020)IINTFP*I090,IINTFP*MEMREQ
 5020        FORMAT(T3,'@GEOPT-F, Not enough core.',/,
     &            T3,'Words required:',I9,' Words available:',I9)
             CALL ERREX
          ENDIF
          CALL PROJPLD(NATOMS,NIRREP,IORDER,ATMASS,Z(I000),
     &                 Z(I010),Z(I020),Z(I030),Z(I040),Z(I050),
     &                 Z(I060))
        ENDIF
C
C CONVERT TO IC AND DUMP TO DISK NOW.
C
        CALL CONVHESS(Z(N3),Z(N4),Z(1),Z(N2),Z(N1),-1)
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N6))
        CALL VIB1(Z(1),Z(N24),Z(N21),Z(N25),Z(N12),Z(N25),Z(N1),Z(N6),
     &            Z(N27), Z(N28), Z(N29), Z(N30),NRX,IStat,Z(N26))
        CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)
        CALL ACES_JA_FIN
        STOP
      ELSEIF(IFLAGS(18).EQ.3
     &       .OR.IFLAGS(18).EQ.4
     &       .OR.IFLAGS(18).EQ.5)THEN
C
C NMR CALCULATION, Note that the last argument is just a
C dummy argument. This is to aviod having two subroutines that
C basically do same things except for one thing. Note the
C above call to PREVIB this arguement has the gradient
C vector. In the following calls Z(N26) does not contain
C any useful info. Ajith Perera 12/2001.
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL NMR1(Z(1),Z(9*NRX+1),Z(10*NRX+1),Z(11*NRX+1),NRX)
        CALL GFNAME('JSO     ',FNAME,ILENGTH)
        INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
        IF(YESNO) CALL JSO1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)
        CALL ACES_JA_FIN
        STOP
      ELSEIF(IFLAGS(18).EQ.8) THEN
C
C NMR COUPLING CONSTANT (SO CONTRIBUTION) CALCULATION
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL JSO1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT',IONE,IJUNK)
        CALL ACES_JA_FIN
        STOP

      ELSEIF(IFLAGS(18).EQ.9) THEN
C
C NMR COUPLING CONSTANT (FC CONTRIBUTION) CALCULATION
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL JFC1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT',IONE,IJUNK)
        CALL ACES_JA_FIN
        STOP
      ELSEIF(IFLAGS(18).EQ.10) THEN
C
C NMR COUPLING CONSTANT (SD CONTRIBUTION) CALCULATION
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL JSD1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT',IONE,IJUNK)
        CALL ACES_JA_FIN
        STOP
      ENDIF
      iflags(18)=itmp
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Non-vibrational calculations really begin here
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Save some information about the xy orientation of the molecule
C     before taking the step.  Used later to give new geometry same
C     orientation.
C
      Err = 0
      IF(OPTRES) CALL GETREC(20,'JOBARC','COORD   ',IINTFP*NX,Q)
      Call SavOri (FPGrp, NAtoms, Q, IAtNum, At1, At2, XYOri, Err)
      IF(OPTRES) GO TO 9876
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Get gradient and Hessian from disk
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     IRecal controls how often we should look for an analytic
C     Hessian.  For a "full-x" optimization, IRecal should be 1.
C
C     NOTE: IRecal is not used elsewhere.  Eventually it should be
C     used to generate input decks that insure the Hessian will be
C     recalculated every IRecal steps.
C
      IHes = -1
      If (IRecal .ne. 0) then
        If (Mod(ncycle, IRecal) .eq. 0 ) IHes = 0
      Else
        If ( ncycle .eq. 0) IHes = 0
      EndIf
C
C If this is the first cycle, and an internal coordinate Hessian
C  is available, then snatch it.
C
      IGRD=0
      IF(IFLAGS2(h_IFLAGS2_geom_opt).EQ.0 .OR. IFLAGS2(138).EQ.1)IGRD=1
      Call READGH (IGRD,IHES,NATOMS,Z(N6),Z(1),Z(N1),Z(N2),IAVGRD,
     &     IAVHES)
      If ( ncycle .eq.0 .and. ihes .eq. 0)THEN
        CALL GETICFCM(Z(N2),NXM6,IICFCM)
        IF(IICFCM.EQ.1)IHES=-1
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Convert gradient to internal coordinates & print it out
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CALL CONVQ(Z(N3), Z(N6), Z(N7), 0, 1)
      CALL FORCEOUT(Z(N7))
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     If we have an analytic Hessian, convert it also...
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      If ( ncycle .eq. 0 .OR. IHES.EQ.0) THEN
        IF(IICFCM.EQ.0)THEN
          CALL CONVHESS(Z(N3),Z(N4),Z(1),Z(N2),Z(N1),1-IAVHES)
        ENDIF
C
C FORCE CURVILINEAR <-- RECTILINEAR TRANSFORMATION IF TS SEARCH
C    AND TRANSFORMATION HAS NOT BEEN TURNED OFF BY EXTERNAL
C    REQUEST (ICURVY>1).
C
        IF((INR.EQ.4.OR.INR.EQ.5).AND.ICURVY.EQ.0)ICURVY=1
        IF(ICURVY.GT.1.AND.INR.EQ.4)WRITE(6,9920)
 9920   FORMAT(T3,' Transformation to curvilinear coordinates has ',
     &       'been turned off.')
        WRITE(6,*) "IAM HERE"
        IF (ICURVY .EQ. 1 .AND. IICFCM .EQ. 0)
     &       CALL TWIDLE(Z(N3),Z(N2),Z(N7),Z(N25),Z(N21),Z(N22),
     &       Z(N23),Z(N24),Z(N26))
        WRITE(6,*) "IAM HERE"
C
C DO NOT WRITE OUT WHOLE HESSIAN UNLESS IT HAS BEEN SPECIFICALLY
C   REQUESTED.  THIS CAN EVENTUALLY BE MADE INTO AN ADJUSTABLE
C   PARAMETER, BUT WILL BE IN CODE FOR NOW.
C
        IF(IPRNT.GE.10) THEN
cSSS             WRITE(6,*) (NOPTI(III), III =1, NXM6)
cSSS             CALL HESSOUT(Z(N2), NXM6, NXM6, 0)
        ENDIF
      ELSE
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     ...otherwise update the Hessian
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     IUPD=2 GIVES BFGS; IUPD=1 GIVES POWELL.
C
        IUPD = 2
        IF (INR .EQ. 4.OR.INR.EQ.5) IUPD = 1
        CALL HUPDATE(Z(N7), Z(N2), Z(1), Z(N9), Z(N10), IUPD)
      END IF
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Compute the step (in internal coordinates
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     SAVE HESSIAN FOR ARCHIVING BEFORE IT IS DESTROYED IN STEP TAKING
C     EMPTY Z(N11) FIRST
C
      CALL ZERO(Z(N10), NXM6 * NXM6)
      CALL ZERO(Z(N1), NXM6)
      CALL VADD(Z(N8), Z(N2), Z(N10), NXM6 * NXM6, 1.D0)
C
C SET FLAG WHICH DETERMINES DIMENSION OF AUGMENTED HESSIAN, WHICH IS
C  USED IN RFA AND TS SEARCHES.  NEED TO BE CAREFUL ABOUT NR AND MANR
C  OPTIONS, HOWEVER, SINCE THESE MAY BE SWITCHED TO RFA IF HESSIAN HAS
C  NEGATIVE EIGENVALUES.  FORCE DIMENSION TO BE (NOPT+1) FOR EVERYTHING
C  BUT TS, WHERE IS MUST BE NOPT.
C
      IAUGHS = 1
      IF(INR .EQ. 4 .OR. INR .EQ. 5) IAUGHS=0
C
      CALL TKSTEP(Z(N7), Z(N2), Z(1), Z(N14), Z(N15), Z(N1), Z(N16),
     &            Z(N17), Z(N18), Z(N19), Z(N20), IAUGHS)
C
      IF ( ncycle .GE. iflags(92)) then
        Write (LuErr,*) '*Maximum number of optimization steps ',
     $       'exceeded.'
        Call ErrEx
      EndIf
 9876 CONTINUE
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Generate new cartesian coordiantes from internals and determine
C the point group.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IF (XYZIN) THEN
         CALL PROCESTEP_XYZ(Z(1), Z(N3), Z(N4), NATOMS, TOTREDNCO)
C
C Here I force, NX=3*NATOMS so that symmetry.F would no complain
C when it try to retrive the Cartesiian coordiantes (the record
C length is 3*NATOMS regardles of the number of internal coordinates.
C Note that NX is set back again to MAX (TOREDUNCO, 3*NATOMS) just
C before archiving. Ajith Perera, 10/2002
C
         NX = 3*NATOMS
      ELSE
         CALL GMETRY
      ENDIF
      CALL SYMMETRY(Z(N11), Z(N12), Z(N13))
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Do some fudging with the orientations
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Eventually, symmetry will always leave the molecule in Cotton's
C     orientation and orientation "fudging" to accomodate specific
C     integral packages will be done here.
C
C     It is possible (in a few points groups like C2v & maybe D2h)
C     for the relative orders of the moments of inertia to switch
C     from one step to the next.  This effectively interchanges
C     the x and y axes, and thus some of the irreps -- wreaking havoc
C     on occupations specified by symmetry irrep.  The solution is
C     to base everything on the orientation of the zeroeth step.
C
c      IF(INR.NE.6)Call FixOri (NAtoms, Q, At1, At2, XYOri)
      Call FixOri (NAtoms, Q, At1, At2, XYOri)
      ITMP=IPRNT
      IF(PGRP.EQ.FPGRP)IPRNT=0
      IOFF=N19+1+120*9+NATOMS*9
      MLTINT=IINTFP*(MEMREQ-N19)
      CALL SYMDRV(Z(IOFF),Z(N19),MLTINT,IOFF,PGRP,'COMP')
      CALL SYMDRV(Z(IOFF),Z(N19),MLTINT,IOFF,FPGRP,'FULL')
      IPRNT=ITMP

C   o print out coordintes and distance matrix
C     Note that all orientation fudging has already been done. These
C     coordinates should match the integral program's except for
C     permutations of equivalent atoms.
      Call GeomOut
      Call ADM

C   o make input decks for the integral program using new coordinates
      If (IntTyp.Eq.0) Then
cYAU         Call MkAbIn (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
cYAU     &                GenBy, Z(N11), IStat)
           write(*,*) '@GEOPT: MkAbin was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
CJDW 9/20/96. Allow for Seward (currently must make MOL file).
      Else If (IntTyp .eq. 1 .or.
     &         IntTyp .eq. 2 .or.
     &         IntTyp .eq. 4) Then
         Call MkVMol (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &                GenBy, Z(N11), IStat, BasNam)
      Else If (IntTyp .eq. 3) Then
cYAU         Call MkCadp (PGrp, NAtoms, NUnique, ZSym, IAtNum, Q,
cYAU     &                GenBy, IStat)
           write(*,*) '@GEOPT: MkCadp was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
      End If
      If (Mod(IStat,2) .eq. 1) then
         Call ErrEx
      Else If (IStat .eq. 8 .OR. IStat .eq. 10) then
         IStat = 0
         Write (LuErr, 9810)
C        Eventually we will want to do more things here...
C        Probably want to change to program to chain to
C        Probably want to write out a CON file or similar
      End IF

c o save the information about this step to disk
c   THE GRADIENT IS IN Z(N7)
c   THE HESSIAN  IS IN Z(N8)
C
      NX = MAX(TOTREDNCO, 3*NATOMS)
      IF (.NOT.OPTRES) CALL ARCHIVE(E,Z(N7),Z(N8),Z(N1),1,Z(N18))
      IF (IPRNT.EQ.999) WRITE(LUOUT,8888) (I,Z(I),I=1,N27)
      CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)
      CALL ACES_JA_FIN

      IF (NCYCLE.GE.IFLAGS(92)) THEN
         WRITE(LUOUT,*)
     &      '@GEOPT: Maximum number of optimization cycles exceeded.'
         WRITE(LUOUT,*)
     &      '        Geometry optimization failed.'
         CALL ERREX
      END IF

c     FORMAT STATEMENTS USED MORE THAN ONCE
 8888 FORMAT(T3,'@GEOPT: Dump of core vector:',/,
     &       T5,'Element',T40,'Value',/,(T5,I6,T36,E15.10))
 9810 FORMAT(T3,'@GEOPT: Standard ACES2 routing impossible. ',/,
     &       T3,'        User must supply gradient and Hessian.')

cYAU - cleans up output
c      STOP
      END

