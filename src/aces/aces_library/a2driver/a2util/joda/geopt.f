
C GEOPT:  A Z-matrix input and geometry optimization module for
C         the ACES2 program system.
C
C By: J.F. Stanton, 7/88
C Modifications for portability by: D.E. Bernholdt 12/88
C Modifications for [name it] by: The ACES Development Team








































































































































































































































































































































































































































































































































      SUBROUTINE GEOPT
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
      LOGICAL IGNORE,OPTRES,YESNO, OPTARC_presnt, FCM_EXIST,
     &        IMPORT_CART_HESS
cYAU : 3/11/03 : VIBRES is only used in this file and for nothing
c      LOGICAL VIBRES
      CHARACTER*(fnamelen) FNAME
C     Allocate more than enough core for the present needs
C     Estimated  by NX-6 => NX, and assume the worst!
C      PARAMETER (MEMREQ = 14*9*MXATMS*MXATMS + 16*3*MXATMS +
C     &           5*9*3*MXATMS + 2*3*MXATMS)
C Audit was carried out and recounted the memory requirements. Also,
C change to MXATMS to MAXREDUNCO so tha both user defined internal
C and redundent intnterals can work (note that No. of REDUNCO 
C = 3*MXATMS, so this is safe for both options), 05/2006, Ajith Perera.
C
      PARAMETER (MEMREQ = 15*MAXREDUNCO*MAXREDUNCO + 
     &           34*MAXREDUNCO) 
      DIMENSION Z(MEMREQ)
C
CJDW  Arrays for the POLYPRT routine.
CJDW  10/23/95. Should not now be needed as we are using POLYPRT0.
C      dimension xyzcor(3*mxatms),grad(3*mxatms),hess(9*mxatms*mxatms)
C      dimension iordco(  mxatms)
C
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
      CHARACTER*7 EXPRT_INT_HESS

cYAU passed down to entry for mkvmol
c     Maximum string length of basis set
      INTEGER BASLEN
      PARAMETER (BASLEN=80)
      character*(baslen) BasNam
      INTEGER TOTREDNCO
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
      LOGICAL XYZIN,PS1EXIST,We_havegeom,Can_do_freq,Anlytic_hessian
     &       ,Geomtry_opt, Do_pes_scan
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      COMMON /INPTYP/ XYZIN,PS1EXIST
      COMMON/RESTART_COM/IGNORE
      COMMON/RESTART2/OPTRES
cYAU      COMMON/RESTART3/VIBRES
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
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C     Use to fix orientation for C2v and D2h from one step to the next
      Integer At1, At2, XYOri
C
      integer genby(mxatms)

c ----------------------------------------------------------------------


      IONE   = 1
      IJUNK  = 1
      XYZTol = 1.0d-6
      SYMTOL = 1.D-6
      DEGTOL = 1.D-10
      ISTAT  = 0
      E      = 0.D0
      IPrnt  = 0
      TOTREDNCO  = 0
      i_havegeom = 0
      i_pes_scan = 0
C
C IF .OPT FILE NOT THERE, READ Z-MATRIX AND GET THINGS SET UP
      call getrec(-1,'JOBARC', 'HAVEGEOM', 1, i_havegeom)
      call getrec(-1,'JOBARC', 'PES_SCAN', 1, i_pes_scan)
      Print*, "The value of HAVEGEOM record:", i_havegeom 
      Print*, "The value of PES_SCAN record:", i_pes_scan
      We_havegeom = .false.
      Can_do_freq = .false.
      Do_pes_scan = .false.
      If (i_havegeom .EQ. 1) We_havegeom = .true.
      If (i_pes_scan .EQ. 1) Do_pes_scan = .true.
      CALL ENTRY(BasNam, We_havegeom, Can_do_freq)

c      inquire(file='ALLDONE',exist=yesno)
c      if( (iflags2(138).eq.2) .and. .not. yesno)iarch=0
C
      call getrec(1,'JOBARC','FNDFDONE',1,iTmp)
      if ((iflags2(138).eq.2).and.(iTmp.eq.0)) iarch=0
      Anlytic_hessian =  (Iflags(54) .EQ. 1)
      Geomtry_opt     =  (Iflags2(5) .NE. 0)
      If (Do_pes_scan) XYZIN = .True.
C
      Print*, "-----After call to Entry-----"
      Print*, "The FNDDONE in Geopt:", iTmp 
      Print*, "The gradient calcs:", iflags2(138)
      Print*, "Hessian calc.,Geo. opt. and iarch:", Anlytic_hessian,
     &         Geomtry_opt, iarch
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
      If (We_havegeom .AND. (Anlytic_hessian .OR. iTmp .EQ. 1)) Then
      Print*, "Reading previous step from JOBARC file"
           Call Getrec(20, 'JOBARC', 'ZMATATMS', 1, NATOMS)
           Call Getrec(20, 'JOBARC', 'LINEAR  ', 1, ILINEAR)
           If (ILINEAR .EQ. 1) Then
               NX   = 3*NATOMS
               NXM6 = NX - 5
           Else
               NX   = 3*NATOMS
               NXM6 = NX -6
           Endif
           IF (iflags(68) .EQ. 1) XYZIN = .TRUE. 
           Call Getrec(20, "JOBARC", "COORD   ", NX*IINTFP, Q)
           Write(6, "(3F10.5)") (Q(I), I=1, NX)
           Call Getrec(-1, 'JOBARC', 'CONCTVTY', NX, NCON)
           Call Getrec(-1, 'JOBARC', 'CORD_INT', NX*IINTFP, R)
CSSS           If (Anlytic_hessian) Call Bohr2angs(R, NX)
           Call Bohr2angs(R, NX) 
           Call Getrec(20, 'JOBARC', 'ATOMMASS', NATOMS*IINTFP,
     &                 ATMASS)
           Call Getrec(20, 'JOBARC', 'ORIENTMT', 9*IINTFP,
     &                 ORIENT)
           Call Getrec(20, 'JOBARC', 'ATOMCHRG', NATOMS, IATNUM)
           Call Getrec(-1, 'JOBARC', 'ICSQUASH', NX, ISQUASH)
           Call Getcrec(20, 'JOBARC', "PTGP    ", 4, FPGRP)
           Call Getcrec(20, 'JOBARC', "ABL_PTGP", 4, BPGRP)
           Call Getcrec(20, 'JOBARC', "CMP_PTGP", 4, PGRP)
           Call Getcrec(20, 'JOBARC', "ZSYM", 5*NATOMS, ZSYM)
           Call Getcrec(-1, 'JOBARC', "INTCNAM", 5*NX, VARNAM)
CSS           IF (.not. Geomtry_opt .And. Iarch .Ne. 1 .And. .not. 
CSS     &          Do_pes_scan) CALL FETCHZ(.TRUE., Z, MAXMEM)
           CALL ZERO(Z, MEMREQ)
           Print*, "The data read from JOBARC"
           Print*, "The NATOMS:", NATOMS
           Write(6,*)
           IF (.not. xyzin) Write(6,*) "Internal coords:"
           IF (.not. xyzin) Write(6,"(3F10.5)") (R(I), I=1, NX)  
           Write(6,*)
           Write(6,*) "The Cartesian coords:"
           Write(6, "(3F10.5)") (Q(I), I=1, NX)
           Write(6,*)
           if (.not. xyzin) Write(6,*) "The connectivities:" 
           if (.not. xyzin) Write(6, "(5I2)")(NCON(I), I=1, NX) 
      Else
c
          Print*, "Entering Fetchz IARCH .NE. 1", IARCH
          IF (iarch .NE.1) CALL FETCHZ(.TRUE., Z, MAXMEM)
      Endif
C
c Marcel Nooijen 2001
c print warning if atoms are connected 3 to 1 and vibrations.
c Abort if resraman calculation
c
      if (iarch .ne. 1) then
c
         ijunk = 0
         call getrec(-1, 'JOBARC', '12SWITCH', ione, ijunk)
         if (ijunk .eq. 1 .and. (iflags2(3) .eq. 1
     $        .and. ivib .ge. 1)) then
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
C The JOBARC record PASS1 is created in symcor ME and can have
C any of the following values: 0 to indicate a first call to
C symcor in a frequency or finite difference geometry optimizations,
C 1 to indicate geometry optimization with numerical gradients
C or frequency calculations with numerical or analytical gradients
C is in the middle of generating required Hessian or gradients.
C Also, 0 is used to indicate that gradients or Hessians 
C generation is complete. When the gradients or Hessian generation
C is complete, the joda Handle the Frequency generation or 
C taking the next geometry step by transfering the execution 
C to the appropriate block (note the comments in the section marked
C as  BRANCHING below). Also, please read comments in SYMCOR ME.
C The second GETREC call returns the length of the record if that
C exists. If the PASS1 record present it is a finite difference 
C calculation. Ajith Perera 03/2003. 
C 
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,IPASS1)
C 
C The following logic is for restarting numerical gradient
C generation. If there is PASS1 record in JOBARC, that means
C symcor has been run at least once. 
C
      CALL GETREC(0,'JOBARC','PASS1',iTmp,IX)
      PS1EXIST=(iTmp.GT.0)
      IF (PS1EXIST.AND.(.NOT.IGNORE)) THEN
cYAU        VIBRES=.true.
         CALL PUTREC(1,'JOBARC','VIBRES',1,1)
      ELSE
cYAU        VIBRES=.false.
      ENDIF
      Write(6,*)
      Print*, "@GEOPT, PASS1,  PS1EXIST=:", IPASS1, PS1EXIST
C
C Allocate memory, Note that in the case of pure Cartesian Optimizations
C the number of internal coordinates is always correspond to
C NX where NX=3*NATOMS. In the case of Redundent internal optimizations
C NX = TOTREDNCO and the fetchz.F dependency, gen_rednintrnls.F
C create the "REDNCORD" record. In the case of user defined internals 
C it is NX-6. Ajith Perera, 09/2002,2008.
C
      CALL GETREC(0, 'JOBARC','REDNCORD', Length, TOTREDNCO)
      IF (iFlags2(5).eq.2) Then 
          NXM6 = NX
      ELSE IF (length .gt.0) THEN
C
C This is equivalent to iFlags2(5).eq.3, but
C the optimizer (ifort on SGI Altix) would not allow me to use 
C it, so the above is the alternative. 
C
          CALL GETREC(-1, 'JOBARC','REDNCORD', 1, TOTREDNCO)
          NXM6 = TOTREDNCO
      ENDIF
C
C The following applies for the user defined internal optimizations.
C Somewhat overaloaded because PS1EXIST control numerical finite
C differencing. Ajith Perera, 08/2008.
C
      IF (iFlags(68).eq.0 .OR. PS1EXIST .or.
     &    iFlags2(5) .eq. 0) NXM6 = MAX(NX-6,1)
C
C FOR DUMMY B-MATRIX (AND CARTESIAN HESSIAN AFTER BUILDB)
C
      N1 = (NX * NX) + 1
C
C FOR B MATRIX (OCCUPIED BY A(TRSPSE) AFTER BUILDB. (++9*MXATMS^2)
C
      N2 = (NX * NX) + N1
C
C FOR G MATRIX (AND INTERNAL COORD HESSIAN AFT(++9*MXATMS^2)
C
      N3 = (NX * NX) + N2
C
C FOR A MATRIX (++9*MXATMS^2)
C
      N4 = (NX * NXM6) + N3
C
C FOR B TRANSPOSE
C Increase the space at N4 so that it can be used to store the Cartesian
C Hessian in get_mopac_hess.F. The N4 location is used as scratch space
C in several other places and the increased space should not have any
C effect on those. Ajith Perera, 06/2005. (++3*MXATMS*(3*MXATMS-6))
C N5 = (NX * NXM6) + N4
C
      N5 = (NX * NX) + N4
C
C FOR SCRATCH ARRAY IN BUILDB (++9*MXATMS^2)
C
      N6 = (NXM6 * 2) + N5
C
C CARTESIAN GRADIENT VECTOR (++2*(3*MXATMS-6))
C
      N7 = NX + N6
C
C INTERNAL COORDINATE GRADIENT VECTOR (++3*MXATMS)
C
      N8 = NXM6 + N7
C
C INTERNAL COORDINATE HESSIAN FOR ARCHIVING (++(3*MXATMS-6))
C
      N9 = (NXM6 * NXM6) + N8
C
C FOR STEP IN UPDATE (++(3*MXATMS-6)^2)
C
      N10 = NXM6 + N9
C
C FOR SCRATCH ARRAYS IN UPDATE (++(3*MXATMS-6))
C
      N11 = ((NXM6 * NXM6) * 3) + N10
C
C FOR SCRATCH ARRAYS IN SYMMETRY (++3*(3*MXATMS-6)^2)
C
      N12 = (3 * NX) + N11
C
      N13 = NX + N12
C
      N14 = NX + N13
C
C FOR CONSTRAINED HESSIAN IN NEWRAP (++5*(3*MXATMS-6))
C
      N15 = (NOPT * NOPT) + N14
C
C FOR CONSTRAINED GRADIENT VECTOR IN NEWRAP (++9*MXATMS^2)
C
      N16 = NOPT + N15
C
C CONSTRAINED ORTHOG. BASIS GRADIENT & HESSIAN IN EFOL (++3*MXATMS)
C
      N17 = NOPT + N16
C
      N18 = NOPT*NOPT + N17
C
C Eigenvector followed on previous iter. (EFOL) (++(9*MXATMS^2+3*MXATMS))
C
      N19 = NOPT + N18
C
C Scratch matrices for EFOL - Used for augmented Hessian in RFA
C and TS searches.  Dimension must be NOPT+1 for RFA.(++3*MXATMS)
C
      N20 = (NOPT+1)*(NOPT+1) + N19
      N21 = (NOPT+1)*(NOPT+1) + N20
      N22 = NX + N21
      N23 = NX + N22
      N24 = NX + N23
      N25 = NX*NX + N24
C
C A Bug fix (Only occured when extended to do Raman intensiies)
C Ajith 08/98, Note the generous allocation. (++(3*9*MXATMS^2 + 3*MXATMS))
C
      N26 = NX*NX + N25
C
C Marcel Need a bigger scracth space. (++9*MXATMS^2)
cSSS N27 = NXM6*NXM6 + N26
C
      N27 = NX*NX + N26
C
C Allocate room for polarizabiliity derivatives Ajith & John, 08/98 (++9*MXATMS^2)
C
      N28 = 9*NX + N27
      N29 = 9*NX + N28
      N30 =   NX + N29
      N31 =   NX + N30
C
C Allocate memory for a scratch array in line search, Ajith, 03/10.
  
      N32 = 6*NX + N31
C 
C The memory for this block (++20*3*MXATMS)
C
      MEMTOP = N32
C
C The total memory required  (9*15*MXATMS^2 +  3*34*MXATMS ) assuming 
C (3*MXATMS-6)=3*MXATMS for user defined internals and 15*MAXREDUNCO^2
C 34*MAXREDUNCO) for redundent internals
C
      IF (MEMTOP .GT. MEMREQ) THEN
        WRITE(LUOUT,5010)MEMTOP,MEMREQ
 5010   FORMAT(T3,'@GEOPT-F, Not enough core.',/,
     &       T3,'Words required:',I9,' Words available:',I9)
        CALL ERREX
        Write(6,*)
      ENDIF
C
C
C Gets information from OPTARC file for subsequent steps (after
C the first call)
C
      IF (iarch .EQ.1) Then
C
C IF archive file (OPTARC) present then immediately branch out to
C the geometry optimization section. OPTARC is made during the very first
C call to joda for optimizations. For analytical gradients optimizations
C iarch is set 1 (after the first joda) while for numerical derivatives 
C it set to 0 until the full set of gardients available when it is 1.
C (Note that this changes from 0 --> 1 during numerical opt. cycles)
C 03/2003, Ajith Perera. 
C
        If (We_havegeom) Then
           NATOMS = NX/3 
           Call Getrec(10, 'JOBARC', 'IFLAGS  ', 100, IFlags)
           IF (iflags(68) .EQ. 1) XYZIN = .TRUE. 
           Call Getrec(20, "JOBARC", "COORD   ", NX*IINTFP, Q)
           Call Getrec(-1, 'JOBARC', 'CONCTVTY', NX, NCON)
           Call Getrec(-1, 'JOBARC', 'CORD_INT', NX*IINTFP, R)
CSSS           Call Angs2bohr(R, NX)
           Call Getrec(20, 'JOBARC', 'ATOMMASS', NATOMS*IINTFP, 
     &                 ATMASS)
           Call Getrec(20, 'JOBARC', 'ORIENTMT', 9*IINTFP, 
     &                 ORIENT)
           Call Getrec(-1, 'JOBARC', 'ICSQUASH', NX, ISQUASH)
           Call Getcrec(20, 'JOBARC', "PTGP    ", 4, FPGRP)
           Call Getcrec(20, 'JOBARC', "ABL_PTGP", 4, BPGRP)
           Call Getcrec(20, 'JOBARC', "CMP_PTGP", 4, PGRP)
           Call Getcrec(-1, 'JOBARC', "INTCNAM", 5*NX, VARNAM)
           CALL ZERO(Z, N11)
           Print*, "The Retrive-2"
           Print*, "Coordiantes:",iflags(68)
           Print*, "The data read from JOBARC"
           Print*, "The NATOMS:", NATOMS
           Print*, "Untested change (undoB2A) on 07/09/06"
           IF (.not. xyzin) Write(6,*) "Internal coords:"
           IF (.not. xyzin) Write(6,"(3F10.5)") (R(I), I=1, NX)
           Write(6,*)
           Write(6,*) "The Cartesian coords:"
           Write(6, "(3F10.5)") (Q(I), I=1, NX)
           Write(6,*)
           If (.not. xyzin) Write(6,*) "The connectivities:"
           If (.not. xyzin) Write(6, "(5I2)")(NCON(I), I=1, NX)
        Else
            CALL RETRIEVE( E, Z(N4), Z(N8), Z(1), Z(N18))
            CALL ZERO(Z, N11)
C
      Write(6,*), "Data read from OPTARC file"
      Write(6,*) 
      Print*, "The COORD COMMON BLOCK/AFTER RETRIVE"
      Write(6,"(3F10.5)"), (Q (I), I= 1, NX)
      Write(6,*)
      Write(6, "(3F10.5)"), (R(I), I= 1, NXM6)
      Write(6,*)
        Endif 
C 
      ELSE
C
C CALL SYMMETRY PACKAGE AND NEW FINDIF SYMMETRY ROUTINES.
C
            Print*, "Before the call to GMETRY, XYZIN", XYZIN
        IF (.NOT. XYZIN) CALL GMETRY(We_havegeom, .TRUE.)
        IF (XYZIN .AND. Iflags(54) .GT. 1 .AND. 
     &      IPASS1 .EQ. 0) CALL GETXYZ
C
C Internal coordinate geo. optimizations, Cartesian Coordinate 
C frequency calcualtions requires following call to SYMMETRY.  
C This call should be avoided for RIC optimizations 
C 03/2003 Ajith Perera. 
C 
        IF (iFlags2(5).lt.3) THEN
           CALL SYMMETRY(Z(N11), Z(N12), Z(N13), .TRUE.)
        ELSE
C
C Allow higher tolerance for Cartesian inputs (perhaps generated from a
C GUI). The parameter DEGTOL is used to determine the degeneracy of the
C inertia tensor which is related to the degeneracy of the point group.
           DEGTOL = 1000*DEGTOL
        END IF
C
CJDW 4/1/97. Write COMPSYOP with maximum possible length at beginning
C of finite difference calculation.
C Also COMPPERM, COMPCLSS. When SYMMETRY=NONE, these two records
C are initialized in symmetry.F.
C
        IF(.NOT.PS1EXIST)THEN
         If (iFlags(60).GT.0) THEN
            CALL ZERO(Z,8*9)
            CALL PUTREC(20,'JOBARC','COMPSYOP',8*9*IINTFP,Z)
            CALL ZERO(Z,8*NATOMS)
            CALL PUTREC(20,'JOBARC','COMPPERM',8*NATOMS,Z)
         Endif
C
            CALL ZERO(Z,8)
            CALL PUTREC(20,'JOBARC','COMPCLSS',8,Z)
         IF (Iflags2(5) .eq. 2 .or. 
     &       Iflags2(5) .ge. 3) THEN
             CALL PUTREC(0,'JOBARC','COORD_OP',IINTFP*3*NATOMS,Z)
         ENDIF
        ENDIF
C
        
        IOFF=1+120*9+NATOMS*9
        MLTINT=IINTFP*(MEMREQ-120)
C
        If (iFlags(60).GT.0) CALL SYMDRV(Z(IOFF),Z,
     &                                 MLTINT,IOFF,FPGRP,'FULL')
        ITMP=IPRNT
        IF(PGRP.EQ.FPGRP)IPRNT=0
        If (iFlags(60).GT.0) CALL SYMDRV(Z(IOFF),Z,
     &                                 MLTINT,IOFF,PGRP,'COMP')
        IPRNT=ITMP
CSSS        CALL GETREC(20,'JOBARC','NOREOCOR',3*NATOMS*IINTFP,
CSSS     &              Q)

        CALL ZERO(Z,MEMREQ)
C
c SB : The if(.not.PS1EXIST) block used to be inside the if(natoms.gt.1)
c block.  I moved it outside so that the call to geomout is made
c even in the case of an atom.
C
        If (NAtoms .gt. 1) THEN
C
C Call the Z-matrix analyzer - requires A matrix from BUILDB.
C XYZIN is true when IFLAGS(68)=1 occur at the same time. So,
C in that sens the following test is a overkill. 03/2003, 
C Ajith Perera.
C 
          IF (.NOT. XYZIN .AND. IFLAGS(68).NE.1) THEN
C
             CALL BUILDB(Z(1),Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),Z(N10))
C
C Analyse only should be called for Geo. optimizations
C
             IF(.NOT.PS1EXIST)Call Analyze (Z(N10), Z(N3))
C
C As of 08/2008 we can also useiFlags2(5).eq.3ne.3)ne.3)
C to test for RIC opts.
C
          ELSE IF (iFlags2(5).ge.3) THEN
C
C Note, this  part is executed only in the first call to joda, and
C no need to re-run since it had alredy been run with the argument set
C to .TRUE, and the following records for the starting  geometry is in
C the JOBARC. Ajith Perera 08/2002.
C
             CALL GETREC(-1, 'JOBARC','REDNCORD', 1, TOTREDNCO)
             LENGTH_GMAT=3*NATOMS*TOTREDNCO
             CALL GETREC(-1, 'JOBARC', 'GMATRIX ',LENGTH_GMAT*IINTFP,
     &                   Z(N3))
             CALL ANALYZE(Z(N10), Z(N3)) 
C
          ELSE IF (iFlags2(5).eq.2 .AND.XYZIN) THEN
C
C This block is for pure Cartesian full optimizations. In this case 
C all the B and G matrices are unit matrices.
C
             CALL SET_2UNIT_MATRIX(Z(N3), 3*NATOMS)
             CALL ANALYZE(Z(N10), Z(N3))     
C
          ENDIF
        endif
C
C Print out the coordinates and distance matrix
C
        IF(.NOT.PS1EXIST)THEN
          Call GeomOut
          if (natoms.gt.1) Call ADM
        EndIf
C
C IF THIS IS A FINITE DIFFERENCE CALCULATION AND ALL POINTS HAVE
C BEEN RUN, THEN CALL VIBRATIONAL ANALYSIS CODE
C Make the appropriate input deck
C
        IF (INTTYP .EQ. 0) THEN
cYAU          Call MkAbIn (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
cYAU     &         GenBy, Z(N11), IStat)
           write(*,*) '@GEOPT: MkAbin was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
CJDW 9/20/96. Allow for Seward (currently must make MOL file).
CADY 5/06/04. Allow for GAMESS (currently must make MOL file).
        ElseIf (IntTyp .eq. 1 .or. IntTyp .eq. 2
     &                        .or. IntTyp .eq. 5
     &                        .or. IntTyp .eq. 4) then
C
C The last argument is added to reformat the basis set to satisfy 
C the requirement for the wave function analysis (see a2proc). The
C reformatting routine read a tmp file, GENB.TMP, created in the
C first MkVmol call and write the ZMAT.BAS file in the new format
C and the MOL file is always created from the the ZMAT.BAS (including 
C the very first run). Tom Watson Jr. and Ajith Perera; 01/2011. 
C
          Call MkVMOL (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &         GenBy, Z(N11), IStat, BasNam)

           if (iflags(45).eq.2) then
              call mknddo(q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &                    GenBy, Z(N11), IStat, BasNam)
           end if
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
        itmp                   = iflags(18)
        iflags(18) = mod(itmp,100)
C
C ****************** BRANCHING Begins **************************
C
C Test to see whether the current calcultion is a single point
C energy, vibrational frequency or geometry optimiztion. Unfortunatley
C some NMR spin-spin coupling stuff (work of Stanton and Gauss not
C Perera's) cramped here (The iflags(h_IFLAGS_props) carry those flags.
C In principle these should have been moved out of joda). If not a single
C point or atomic, create archive file OPTARC, else print a message and
C move on to next ME.
C 
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
c
c inquire(file='ALLDONE',exist=yesno)
c if(.not.yesno)then, if(.not. PS1EXIST .or. PS1EXIST .and. IPASS1 .eq. 1)
C These lines were obsoleted on 03/2003 Ajith Perera
C
C The logic here is important. The variable IPASS1 can be 0 and 1
C (see the comments under PASS1 JOBARC record above). 
C Geometry Optimizations with analytical gradients reach this far
C only in the first call to joda and the condition is true. In 
C subsequent steps, it will see the presence of OPTARC file and
C branch out to the geometry optimization section (below). It can
C exit via SUMMARY.F when the convergence is achieved. 
C 
C For geometry optimizations with numerical gradients the conditional
C is true until the full set of gradients are created (IPASS1 is 
C 1 until the full set of gradients are created then it is set to 0
C in symcor ME). When the full set of gradients are there, iarch get
C set to 1 (OPTARC is present), and the execution get pass over to 
C the geometry optimization logic. It can exit via SUMMARY.F. 
C
C For frequency calculations, untill the full Hessian is available the
C the conditional is true (IPASS1 is 1 PS1EXIST is true). When the 
C conditional is false, execution is transfered to VIB1.F. The 
C program exit is from the main program GEOPT.F. Ajith Perera 03/2003. 
C 
C Note that during the first run IPASS1 is 0 and PS1EXIST is false. 
C 
        iflags(18)=itmp
        if(IPASS1 .eq. 1 .or. .not. PS1EXIST)then
          IF(IPRNT.EQ.999)WRITE(LUOUT,8888)(I,Z(I),I=1,N27)
          CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)
          return
        endif
      END IF
 9820 Format (' @GEOPT-W, Archive file not created for ',A,
     $     ' calculation.')
C 
C ************************ BRANCHING ENDS *********************************
C
C CALL ROUTINE TO BUILD B MATRIX, Only for geometry optimizations
C vibrational Frequency calculations. Also, here, XYZIN and 
C IFLAGS(68) convey the same information. 03/2003, Ajith Perera
C
      IF(.NOT. XYZIN .AND. IFLAGS(68).NE.1) THEN
C
        CALL BUILDB(Z(1), Z(N1), Z(N2), Z(N3), Z(N4), Z(N5), Z(N10)) 
C
      ELSE IF (TOTREDNCO.NE.0) THEN
C
C This block is executed for each iteration during a optimiation and
C requires a call to FETCHZ with argument set to .FALSE. to create
C B and A matrices for the new coordinates using the IRC's genereated
C in the first call to FETCHZ. Ajith Perera 08/2002
C
        CALL FETCHZ(.FALSE., Z, MAXMEM)
        CALL GETREC(-1, 'JOBARC','REDNCORD', 1, TOTREDNCO)
        LENGTH_GMAT = 3*NATOMS*TOTREDNCO
        CALL GETREC(-1, 'JOBARC', 'GMATRIX ',LENGTH_GMAT*IINTFP,
     &              Z(N3))
        CALL GETREC(-1, 'JOBARC', 'BMATRIXT', LENGTH_GMAT*IINTFP,
     &              Z(N4))
C
       ELSE IF (iFlags2(5).eq.2 .AND.XYZIN) THEN
C
C This block is for pure Cartesian full optimizations. In this case
C all the B and G matrices are unit matrices.
C
             CALL SET_2UNIT_MATRIX(Z(N4), 3*NATOMS)
             CALL PUTREC(20,'JOBARC','BMATRIX ',9*NATOMS*NATOMS*
     &                   IINTFP, Z(N4))
             CALL SET_2UNIT_MATRIX(Z(N3), 3*NATOMS)
             CALL PUTREC(20, 'JOBARC', 'AMATRIX', 9*NATOMS*NATOMS*
     &                   IINTFP, Z(N3))

      ENDIF
C
C************************* Frequency/Rate Cons/NMR starts **********
C 
C This start the vibrational Frequency,  Stanton Gauss NMR stuff
C and Rozeann Steckler's rate constant stuff. 
C
      itmp                   = iflags(18)
      iflags(18) = mod(itmp,100)
C
      IF (Can_do_freq .OR.IFLAGS2(3).NE.0)THEN
C
CJDW 1/ 5/95 Write polyrate data for Rozeanne Steckler. This almost
C certainly will not work if there are dummy atoms present.
C In addition, it should only be called when we do analytical
C hessians, as in the finite difference case it will not get
C the coordinates, energy, and gradient for the reference
C geometry. Also, POLYRATE must be ON.
C
CJDW 10/23/95 Change from POLYPRT to POLYPRT0 in standard code.
C
        IF(IFLAGS(54).EQ.1.AND.
     &     IFLAGS2(113).EQ.1)THEN
          CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                  .TRUE.,.TRUE.,.TRUE.,.TRUE.)
        ENDIF
C
        IHES = 0
        IF(IVIB.EQ.3.OR.IFLAGS2(3).NE.0) IHES = 1
C
        Call READGH (0,IHES,NATOMS,Z(N6),Z(1),Z(N1),Z(N2),IAVGRD,
     &       IAVHES)

       Print*, "The exact Hessian after reading in READGH."
       CALL OUTPUT(Z(1), 1, 3*NATOMS, 1, 3*NATOMS, 3*NATOMS,
     &             3*NATOMS, 1) 
C
C
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
             WRITE(LUOUT,5020)IINTFP*I070,IINTFP*MEMREQ
 5020        FORMAT(T3,'@GEOPT-F, Not enough core.',/,
     &            T3,'Words required:',I9,' Words available:',I9)
             CALL ERREX
          ENDIF
          CALL PROJPLD(NATOMS,NIRREP,IORDER,ATMASS,Z(I000),
     &                 Z(I010),Z(I020),Z(I030),Z(I040),Z(I050),
     &                 Z(I060))
        ENDIF
C
C CONVERT TO IC AND DUMP TO DISK NOW. This block is executed
C for VIB=EXACT, FINDIF. 
C When the last argument is set to -1, the CONVHESS routine writes out the
C FCMINT file (note that any other value for this argument will not work).
C The character of the FCMINT (which coordinates it is on) depends on
C the input coordinate sytem. Note that:
C
C  1. The redundant internal coordinates (RIC) input - The incoming Hessian
C     is in Cartesian and the outgoing Hessian is in RIC
C     (not implemented yet). So, the outgoing Hessian is in Cartesians.
C
C  2. The user defined internals (ZMAT) input - The incoming Hessian is
C     in Cartesian and the outgoing Hessian is in user defined
C     internals.
C
C  3. Cartesian input - The incoming Hessian is in Cartesian
C     and the outgoing Hessian is in Cartesian.
C
C A cautionary note about the transformation of Cartesian Hessian
C to internal Hessian. Inside CONVHESS it is assumed that the
C curvature of the transformation of internal-cartesian or vice versa
C is zero. Notice the call to TWIDLE below when this Hessian
C is used as an initial guess for a transition state searches where this
C condition is no longer valid. In TWIDLE the terms that arise from
C the curvature of the transformation is evaluated and added to the
C internal Hessian. Ajith Perera, 06/2004.
C
        CALL CONVHESS(Z(N3),Z(N4),Z(1),Z(N2),Z(N1),-1)
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N6))
        CALL VIB1(Z(1),Z(N24),Z(N21),Z(N25),Z(N12),Z(N25),Z(N1),Z(N6),
     &            Z(N27), Z(N28), Z(N29), Z(N30),NRX,IStat,Z(N26))
        CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)
c YAU: leave the Hessian so we can re-run xjoda with ISOMASS files
c (A lot of single-point stuff will break, but those calcs made no
c sense without reference info anyway.)
c        CALL ACES_JA_TRUNCATE('NEXTGEOM',1)
        call rm_backup
        return
C
C Vibrational Frequency Calculations Ends!!!
C
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
        return
      ELSEIF(IFLAGS(18).EQ.8) THEN
C
C NMR COUPLING CONSTANT (SO CONTRIBUTION) CALCULATION
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL JSO1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT',IONE,IJUNK)
        return

      ELSEIF(IFLAGS(18).EQ.9) THEN
C
C NMR COUPLING CONSTANT (FC CONTRIBUTION) CALCULATION
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL JFC1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT',IONE,IJUNK)
        return
      ELSEIF(IFLAGS(18).EQ.10) THEN
C
C NMR COUPLING CONSTANT (SD CONTRIBUTION) CALCULATION
C
        CALL PREVIB(NRX,Z(1),Z(N11),Z(N26))
        NRX=NRX/3
        CALL JSD1(Z(1),Z(NRX*9+1),NRX)
        CALL PUTREC(20,'JOBARC','JODAOUT',IONE,IJUNK)
        return
C
C At the moment we don't need anything here. All we did 
C was to create a numerical energy gradient. What you
C want to do with is up to the person who wants it. 
C
      ELSEIF((iflags2(138).eq.2) .AND.
     &       (iFlags2(5) .eq.0)) THEN 
      
           RETURN
      ENDIF
C
C*************** Frequency/Rate Cons/NMR Ends ************************
C*************** Non-vibrational calculations really begin here ******
C
C Save some information about the xy orientation of the molecule
C before taking the step.  Used later to give new geometry same
C orientation.
C
      Err        = 0
      iflags(18) = itmp
      IF(OPTRES) CALL GETREC(20,'JOBARC','COORD   ',IINTFP*NX,Q)
      Call SavOri (FPGrp, NAtoms, Q, IAtNum, At1, At2, XYOri, Err)
      IF(OPTRES) GO TO 9876
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Get gradient and Hessian from disk
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C IRecal controls how often we should look for an analytic
C Hessian.  For a "full-x" optimization, IRecal should be 1.
C IRecal get initilized in mkvmol.F (??). 
C
C NOTE: IRecal is not used elsewhere.  Eventually it should be
C used to generate input decks that insure the Hessian will be
C recalculated every IRecal steps.
C
C With the introduction of the INIT_HESSIAN keyword, we are now
C able to control the make up of the initial guess Hessian.
C Follow the notes below to see how this new key token is integrated
C to the structure that existed before without breaking it.
C Naturally, that introduces redundancies and since they
C are harmless, I will leave them as they are. Ajith Perera, 02/2005.
C
      INT_HESS = IFLAGS2(8)
      IF      (INT_HESS.EQ.0 .AND. IRECAL .LT. 0) THEN
                                   EXPRT_INT_HESS = "SPECIAL"
      ELSE IF (INT_HESS.EQ.1) THEN
                                   EXPRT_INT_HESS = "FCMINT "
      ELSE IF (INT_HESS.EQ.2) THEN
                                   EXPRT_INT_HESS = "MOPAC  "
      ELSE IF (INT_HESS.EQ.3) THEN
                                   EXPRT_INT_HESS = "EXACT  "
      END IF
C
      Write(6,*)
      Print*, "The COORD COMMON BLOCK/at opt. start:"
C      Write(6,*)
C      Write(6,"(3F10.5)"), (Q(I), I= 1, 3*NATOMS)
      Write(6,*)
      Write(6,"(3F10.5)"), (R(I), I= 1, NXM6)
      Write(6,*)
C
C This is for those who don't want to give explicit instructions. 
C By specifying IRECAL=EVAL_HESS=n, they are telling the program
C to recalculate the SCF Hessian every n cycles, but failed to 
C specify INIT_HESSIAN=EXACT, lets do it for them here. 
C When n=0 SCF hessian generated only once at the begining 
C Ajith Perera, 07/2006. 
C
      IF (IRECAL .GE. 0 .AND. INT_HESS .EQ. 0) 
     &    EXPRT_INT_HESS = "EXACT  "  
C       
C Some assertions; Anthony is going to love this.
C
      If (EXPRT_INT_HESS .EQ. "FCMINT ") THEN
          CALL GFNAME('FCMINT  ',FNAME,ILENGTH)
          INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=FCM_EXIST)
          IF (.NOT. FCM_EXIST) THEN
             Write(6, 01)
  01         Format('INIT_HESS=FCMINT requiers a startating Hessian:',
     &              'starting Hessian matrix (FCMINT) does not exist')
             CALL ERREX
           ENDIF
      ENDIF
C
      IF (EXPRT_INT_HESS .EQ. "SPECIAL" .AND.  IRECAL .GT. 0) THEN
         Write(6, 02)
  02     Format("Inconsistent input parameters: INIT_HESS=SPECIAL ", 
     &          "and EVAL_HESS>0 is not allowed!")
         CALL ERREX
      ENDIF
C
      IHes = -1
      If (IRecal .gt. 0) then
        If (Mod(ncycle, IRecal) .eq. 0 ) IHes = 0
      Else If (IRecal .eq. 0) then
        If (ncycle .eq. 0) Then 
            IHes = 0
          IICFCM = 0
        Endif
      
      EndIf
C
C If this is the first cycle, and an internal coordinate Hessian
C is available, then snatch it. This comment is false. If there
C is an exact Cartesian Hessian then take it. 08/2008.
C
      IGRD = 0
      IF(IFLAGS2(5).EQ.0 .OR.
     &   IFLAGS2(138).EQ.2)IGRD=1
C
      Call READGH (IGRD,IHES,NATOMS,Z(N6),Z(1),Z(N1),Z(N2),IAVGRD,
     &     IAVHES)
C
      Print*, "The Cart. Hessian in after reading from READGH"   
      Print*, "IAVHES flags: Cartesian exact Hess. read",IAVHES,IHESS
      Print*, "If the Hessian is from VIB=EXACT, the trans/rot" 
      Print*, "contaminants are included. If VIB=FINDIF it is"
      Print*, "projected."
      CALL OUTPUT(Z(1), 1, 3*NATOMS, 1, 3*NATOMS, 3*NATOMS,
     &            3*NATOMS, 1)
C
C The following comments deal with what is going on in GETICFCM.
C
C During the first cycle (ncycle=0, ihes=0), read the Internal
C coordinate Hessian from FCMINT file, if FCMINT is avilable,
C then set the IICFCM variable to 1. With the introduction of
C the INIT_HESSIAN keyword, we now know beforehand whether we are
C reading the Hessian from FCMINT file. The function IICFCM becomes
C redundant with EXPRT_INT_HESS = "FCMINT", but I will keep it intact
C since that redundancy does no harm.
C The FCMINT is generated from a prior vibrational calculation.
C Note that the Hessian read in GETICFCM is in Internal coordinates.
C
C IAVHES is 1 if there was an analytic Hessian and zero otherwise.
C This is equivalent to EXPRT_INT_HESS = "EXACT" and once again
C by keeping both does no harm. If there is an analytic Hessian then
C that Hessian is in Cartesian.
C
C Note that the logic here forces us to read FCMINT file regardless
C of whether analytic Hessian has already been read (to allow the
C user to have control over the starting Hessian).
C Ajith Perera, 07/2003.
C
C Let's assume that we have the exact Hessian and
C it is already read by READGH.F, but for some unknown
C reason, a user wanted his/her own Hessian instead. As long
C as the user put the FCMINT file in the workdir, he/she
c can accomplish that goal. The EXPRT_INT_HESS "EXACT" and
C the "FCMINT" are the redundancies.
C
C This is a messy code, in order to have a Cartesian Hessian if
C XYZIN is on, which we are writing to FCMINT (which is also
C used for internal (ZMAT) Hessian) during a vibrational
C frequency calculation with Cartesian input, is read.
C If the input is in cartesian, the Hessian written is in cartesian
C and the variables are set up to accomodate that.
C Luis Galiano, 07/2003. (I have a hard time understanding this
C comment, but I will leave it. Ajith Perera, 02/2005.)

      IICFCM=0
C
C.and. ihes .eq. 0 .and. is removed. 06/08.
C
      If (ncycle .eq.0 .and. (EXPRT_INT_HESS .eq. "FCMINT ")) THEN
C
         IF (.NOT.XYZIN) THEN
C
C Read internal Hessian from FCMINT file!
C
            CALL GETICFCM(Z(N2),NXM6,IICFCM)
            IF(IICFCM.EQ.1)IHES=-1
C
         ELSE
C
C Read Cartesian Hessian from FCMINT file!
C
            If (iFlags2(5).eq.2) Then
                CALL GETICFCM(Z(N2),3*NATOMS,IICFCM)
            Else
               CALL GETICFCM(Z(1),3*NATOMS,IICFCM)
            Endif
            IF(IICFCM.EQ.1)IHES=-1
C
         END IF
      EndIf
      Print*, "The Cartesian Hessian after reading from GETICFCM"
      if (ncycle .eq.0 .and. (EXPRT_INT_HESS .eq. "FCMINT ")) then
      if (iFlags2(5).eq.2) Then
      CALL output(Z(N2), 1, 3*NATOMS, 1, 3*NATOMS,  3*NATOMS, 
     &            3*NATOMS, 1)
      else
      CALL output(Z(1), 1, 3*NATOMS, 1, 3*NATOMS,  3*NATOMS,
     &            3*NATOMS, 1)
      endif
      endif
      Write(6,*)
      Print*, "The COORD COMMON BLOCK/at opt. start:"
C      Write(6,*)
C      Write(6,"(3F10.5)"), (Q(I), I= 1, 3*NATOMS)
      Write(6,*)
      Write(6,"(3F10.5)"), (R(I), I= 1, NXM6)
      Write(6,*)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Convert gradient to internal coordinates & print it out
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      IF (iFlags2(5).eq.2) Then
          CALL DCOPY(NOPT, Z(N6), 1, Z(N7), 1)
      ELSE
         CALL CONVQ(Z(N3), Z(N6), Z(N7), 0, 1)
      ENDIF
      CALL FORCEOUT(Z(N7))
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C If it is the first cycle and there is no Hessian read so far (IAVHES
C is equal to zero, nothing was retrieved from the call to READGH)
C and if a Hessian is needed (IHES=0), then create an approximate
C internal coordinate Hessian (1.0 in the diagonal and 0.0 in
C the off-diagonal). Once again the control is in the final argument
C which takes up the value 1 (1-IAVHESS see below). This is the
C situation when EXPRT_INT_HESS="SPECIAL". The special guess is made
C inside CONVHESS.F.
C
C If analytic Cartesian Hessian was read in READGH (IAVHESS=1), convert
C that to internal coordinates. The conditional based on IICFCM is
C very important. IICFCM is zero if no FCMINT file was present and
C so the Hessian is not read from the FCMINT file. If it is read from
C the FCMINT, it is in Internal Coordinates and the first stage of the
C transformation to Cartesian is not required and hence the call to
C CONVHESS should be avoided.
C  Ajith Perera 07/2003.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      If ( ncycle .eq. 0 .OR. IHES.EQ.0) THEN
        IF(IICFCM.EQ.0)THEN

C IAVHES = 1 (the value of the final argument to CONVHESS is Zero).
C If the Hessian is analytic (EXPRT_INT_HESS = "EXACT"), then it is in
C Cartesian. CONVHESS does the first stage of the transformation to
C internals, i.e. the part that does not involve the gradient of the B matrix.
C
C For EXPRT_INT_HESS = "MOPAC", the NDDO Hessian is read from the
C HINP file and transformed to internal coordinates.
C
C If no Hessian is read (IAVHES=0 or EXPRT_INT_HESS = "SPECIAL")
C and the final argument in CONVHESS is 1. Then CONVHESS generates
C a very crude starting Hessian (this is when no Hessian is read
C either from FCMINT or through READGH and the "MOPAC" Hessian is
C not available). Then a call to TWIDLE is made if the calculation is a
C transition state search. I can identify what goes in TWIDLE
C as the second part of the transformation of Cartesian Hessian
C to internal Hessian (the one that involves the derivative of B-matrix).
C But then, should it not be called for the EXPRT_INT_HESS = "EXACT",
C EXPRT_INT_HESS = "MOPAC", and EXPRT_INT_HESS = "FCMINT"? (All these
C are Cartesian Hessian to start with and CONVHESS only does the
C first stage of the transformation and I thought TWIDLE is
C the one that does the second stage - that is wrong). Whatever
C happens in TWIDLE only effects the Internal Hessian generated
C by an educated guess. In my view, in all three other options (EXACT, MOPAC,
C and FCMINT), the transformation to internals assumes that the
C gradient at that point is zero. This might be fair, since a harmonic
C frequency normally assumes that anyway. Ajith Perera, 02/2005.
C
C If the FCMINT file is generated via a vibrational frequency calculation,
C there is no reason to be here (EXPRT_INT_HESS = "FCMINT" since during
C the frequency calculation, transformation to internals is carried
C out by controlling the final argument to CONVHESS (-1). See above
C vibrational frequency block for further explanation. The overloading
C CONVHESS (generating approximate Hessian in internals, first stage of
C the transformation to internal, and writing the FCMINT file to disk)
C has caused enormous amounts of pain when tried to change or
C understand what is going on there. OVERLOADING should be avoided at
C ALL COST. Ajith Perera 05/2004.

           IF (EXPRT_INT_HESS.EQ."MOPAC  ".AND.IAVHES.EQ.0) THEN
C
C MOPAC Hessian in Cartesians Z(1) or internals Z(N2) if no exact Hessian
C is read!
C 
              CALL GETREC(1,'JOBARC','NREALATM',IONE,NREAL)
              CALL GETREC(1,'JOBARC','LINEAR  ',IONE,ILINEAR)
              NREAL3   = 3*NREAL
              NREAL3M6 = 3*NREAL - 6 + ILINEAR
              CALL GET_MOPAC_HESS(Z(N3), Z(N4), Z(N5), Z(1), Z(N2),
     &                            Z(N1), NATOMS, NX, NXM6, NREAL,
     &                            NREAL3, NREAL3M6)
       Write(6,*)
       Print*, "The Cart. Hessian generated by MOPAC"
       CALL OUTPUT(Z(1), 1, 3*NATOMS, 1, 3*NATOMS, 3*NATOMS,
     &             3*NATOMS, 1)
C
           ELSE IF (EXPRT_INT_HESS.EQ."SPECIAL".AND.IAVHES.EQ.0) THEN
C 
C Educated guess Hessian in internal or redundent internals Z(N2)!
C
       Write(6,*)
       Print*, "Entering convhess:form an educated guess4hess"
       Write(6,*)
              CALL CONVHESS(Z(N3),Z(N4),Z(1),Z(N2),Z(N1),1-IAVHES)
C
           ELSE IF (.NOT. XYZIN .AND. IAVHES.EQ.1) THEN
C
C Cartesian exact Hessian is read Z(1), do part of the transformation 
C that does not depends on the derivative of the B matrix!
C
       Print*, "Entering convhess:piece of cart2int transformation"
       Print*, "IAVHES: ", IAVHES
                 cALL CONVHESS(Z(N3),Z(N4),Z(1),Z(N2),Z(N1),-IAVHES)
C
           END IF

C       ENDIF(IICFCM.EQ.0)
        ENDIF
C
C If redundant internals are used and Cartesian Hessian exists, then
C convert it into redundant internals. This is only needed for
C transition state searches. Let's document the conditional to call
C CHESS:
C
C   TS_SEARCH - simple, applies only for transition state searches.
C   XYZIN .AND. IAVHES.EQ.1   -- RIC (XYZIN=TRUE) and analytic Cartesian
C                                Hessian is available (EXACT).
C   XYZIN .AND. IICFCM .EQ. 1 -- Hessian is read from FCMINT. Note
C                                that the approx. Hessian is generated
C                                from a frequency calculation and it
C                                is in Cartesian.
C   IICFCM .EQ.0              -- An educated guess at the initial Hessian
C                                in internals is made in CONVHESS.F
C
C What is happening in TWIDLE should have been done somewhere else.
C It is actually the part that involves gradient and derivative of B
C matrix in the transformation of Cartesian Hessian to Internal
C Hessians. This comment may be wrong. The TWIDLE is not called
C for when the initial Hessian is EXACT, MOPAC, or read from the
C FCMINT file. If the comment is correct, that should have been
C the case. The call to TWIDLE is made for educated guess Hessian
C (SPECIAL) which is in internal to begin with. Ajith Perera, 02/2005.
C
C In the case of RIC TS searches, this part of the transformation
C (Cartesian to internal) is taken care of in CART2INT_HESS).
C
        IMPORT_CART_HESS = (EXPRT_INT_HESS .EQ. "FCMINT " .OR.
     &                      EXPRT_INT_HESS .EQ. "MOPAC  " .OR. 
     &                      EXPRT_INT_HESS .EQ. "EXACT  " ) 
CSSS        Write(6,*) "IMPORT_CART_HESS", EXPRT_INT_HESS 
C
        IF (IMPORT_CART_HESS .AND. XYZIN .AND.
     &      iFlags2(5).ge.3) THEN
       Print*, "The RIC: Importing Cartesian Hessian",
     &          IMPORT_CART_HESS
       Print*, "The Cart. Hessian in geopt-before CART2INT"
       CALL OUTPUT(Z(1), 1, 3*NATOMS, 1, 3*NATOMS, 3*NATOMS,
     &             3*NATOMS, 1)
C
           CALL CART2INT_HESS(Z(N7),Z(1),TOTREDNCO,NATOMS,Z(N2))
        END IF
C
C FORCE CURVILINEAR <-- RECTILINEAR TRANSFORMATION IF TS SEARCH
C AND TRANSFORMATION HAS NOT BEEN TURNED OFF BY EXTERNAL
C REQUEST (ICURVY>1).
C
CSSS        IF (.NOT.XYZIN) THEN
C
        IF (iFlags2(5).eq.1) Then
           IF (IMPORT_CART_HESS .AND. ICURVY.EQ.0) ICURVY = 1 
           IF (IMPORT_CART_HESS .AND. ICURVY .GT. 1) WRITE(6,9920)
 9920      FORMAT(T3,' Transformation to curvilinear coordinates ',
     &               'has been turned off.')
C
       Print*, "Importing Cartesian Hessian",IMPORT_CART_HESS,
     &          ICURVY
       Print*, "The Partialy trans. Cart. Hessian in geopt-before" 
       Print*, "TWIDLE"
       CALL OUTPUT(Z(N2), 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
           IF (ICURVY .EQ. 1 .AND. IMPORT_CART_HESS) THEN
              CALL TWIDLE(Z(N3),Z(N2),Z(N7),Z(N25),Z(N21),
     &                    Z(N22),Z(N23),Z(N24),Z(N26))
           END IF
        END IF
C
C DO NOT WRITE OUT WHOLE HESSIAN UNLESS IT HAS BEEN SPECIFICALLY
C REQUESTED.  THIS CAN EVENTUALLY BE MADE INTO AN ADJUSTABLE
C PARAMETER, BUT WILL BE IN CODE FOR NOW.
C
CSSS        IF(IPRNT.GE.10) THEN
          CALL HESSOUT(Z(N2),NXM6,NXM6,0)
CSSS        ENDIF
C
C The logic for cycles other than the first cycle begins here.
C Also, note that regeneration of a new Hessian rather than an
C update is requested by the user once every N cycles. During
C the regeneration cycles, it follows the same path as the first
C cycle.
C
      ELSE
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C otherwise update the Hessian, flags(18)=itmpIUPD=2 GIVES BFGS; 
C IUPD=1 GIVES POWELL.
C The Hessian update routine was modified to
C improve the readability and to add new methods. The new additions
C include Murtagh-Sargent update and the mixture, Phi*Powel + (1 - Phi)
C Murtagh-Sargent  (JCC vol 15 pg 1 1994). For minima searches, the
C default is set to BFGS update, and for transition state searches, the
C mixture of Powel and Murtagh-Sargent using Bofill's weighted update
C (see above formula). Ajith Perera, 11, 2004.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
        IUPD = IFLAGS2(7)
        CALL HUPDATE(Z(N7), Z(N2), Z(1), Z(N9), Z(N10), IUPD)
      END IF
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Compute the step (in internal coordinates
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C SAVE HESSIAN FOR ARCHIVING BEFORE IT IS DESTROYED IN STEP TAKING
C EMPTY Z(N11) FIRST
C
      CALL ZERO(Z(N10), NXM6 * NXM6)
      CALL ZERO(Z(N1), NXM6)
      CALL ZERO(Z(N8), NXM6*NXM6)
      CALL VADD(Z(N8), Z(N2), Z(N10), NXM6 * NXM6, 1.D0)
C
C SET FLAG WHICH DETERMINES DIMENSION OF AUGMENTED HESSIAN, WHICH IS
C USED IN RFA AND TS SEARCHES.  NEED TO BE CAREFUL ABOUT NR AND MANR
C OPTIONS, HOWEVER, SINCE THESE MAY BE SWITCHED TO RFA IF HESSIAN HAS
C NEGATIVE EIGENVALUES.  FORCE DIMENSION TO BE (NOPT+1) FOR EVERYTHING
C BUT TS, WHERE IS MUST BE NOPT.
C
      IAUGHS = 1
CSSS      IF(INR .EQ. 4 .OR. INR .EQ. 5) IAUGHS=0
C
C Copy the updated Hessian for approximate frequencies in the final
C iteration. We need to keep a copy of the Hessian because tkstep will
C symmetrize it. Ajith Perera. 05/2005.
C
      CALL DCOPY(NXM6*NXM6, Z(N2), 1, Z(N25), 1)
      CALL TKSTEP(Z(N7), Z(N2), Z(1), Z(N14), Z(N15), Z(N1), Z(N16),
     &            Z(N17), Z(N18), Z(N19), Z(N20), Z(N3), Z(N25),
     &            Z(N4), Z(N31), IAUGHS)
C
C Save the energy of the current step for the trust
C radius step size control algorithm. 11/04, Ajith Perera.
      call getrec(1,'JOBARC','TOTENERG',IINTFP,TotEng)
      call putrec(1,'JOBARC','OLDENERG',IINTFP,TotENg)

C Save the gradient of current step and the geometry for the
C line search algorithm. 12/04, Ajith Perera
c      call dscal(nopt, -1.0D0, Z(nopt+1), 1)
c      call daxpy(nopt, 1.0D0, Z(1),  1, Z(nopt+1), 1)
      call putrec(1,'JOBARC','OLDGRADS',NXM6*IINTFP,Z(N7))
c      call putrec(1,'JOBARC','OLDGEOMT',NOPT*IINTFP,Z(NOPT+1))

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
C
      IF (XYZIN .AND. (iFlags2(5).ge.3)) THEN
CSSS         LENGTH_GMAT=3*NATOMS*TOTREDNCO 
CSSS         CALL GETREC(20, 'JOBARC', 'BTGMIN ',LENGTH_GMAT*IINTFP,
CSSS     &               Z(N3))
         CALL PROCESTEP_XYZ(Z(1), Z(N3), TOTREDNCO)
C
C Here I force, NX=3*NATOMS so that symmetry.F would no complain
C when it try to retrive the Cartesiian coordiantes (the record
C length is 3*NATOMS regardles of the number of internal coordinates).
C Note that NX is set back again to MAX(TOREDUNCO, 3*NATOMS) just
C before archiving. Ajith Perera, 10/2002
C
         NX = 3*NATOMS
C
      ELSE IF ((iFlags2(5).eq.2)) THEN
C
         CALL DCOPY(NX, R, 1, Q, 1) 
         CALL PUTREC(20, 'JOBARC', 'COORD_OP', IINTFP*NX, Q)
C
           Write(6,*)
           Write(6,"(a,a)"), "Updating the new 3*natoms geometries",
     &     " for Cartesian optimizations"
           Write(6, "(a)"),"The Cartesian coords:"
           Write(6, "(3F10.5)"), (Q(I), I=1, NX)
C
      ELSE 
         CALL GMETRY(We_havegeom, .TRUE.)
      ENDIF
      CALL SYMMETRY(Z(N11), Z(N12), Z(N13), .TRUE.)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Do some fudging with the orientations
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Eventually, symmetry will always leave the molecule in Cotton's
C orientation and orientation "fudging" to accomodate specific
C integral packages will be done here.
C
C It is possible (in a few points groups like C2v & maybe D2h)
C for the relative orders of the moments of inertia to switch
C from one step to the next.  This effectively interchanges
C the x and y axes, and thus some of the irreps -- wreaking havoc
C on occupations specified by symmetry irrep.  The solution is
C to base everything on the orientation of the zeroeth step.
C
c      IF(INR.NE.6)Call FixOri (NAtoms, Q, At1, At2, XYOri)
C
      Call FixOri (NAtoms, Q, At1, At2, XYOri)
      ITMP=IPRNT
      IF(PGRP.EQ.FPGRP)IPRNT=0
      IOFF=N19+1+120*9+NATOMS*9
      MLTINT=IINTFP*(MEMREQ-N19)
      If (iFlags(60).GT.0) CALL SYMDRV(Z(IOFF),Z(N19),
     &                               MLTINT,IOFF,PGRP,'COMP')
      If (iFlags(60).GT.0) CALL SYMDRV(Z(IOFF),Z(N19),
     &                               MLTINT,IOFF,FPGRP,'FULL')
      IPRNT=ITMP
C
C print out coordintes and distance matrix
C Note that all orientation fudging has already been done. These
C coordinates should match the integral program's except for
C permutations of equivalent atoms.
C
      Call GeomOut
      Call ADM

C  make input decks for the integral program using new coordinates
C
      If (IntTyp.Eq.0) Then
cYAU         Call MkAbIn (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
cYAU     &                GenBy, Z(N11), IStat)
           write(*,*) '@GEOPT: MkAbin was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
CJDW 9/20/96. Allow for Seward (currently must make MOL file).
CADY 5/06/04. Allow for GAMESS (currently must make MOL file).
      Else If (IntTyp .eq. 1 .or.
     &         IntTyp .eq. 2 .or.
     &         IntTyp .eq. 5 .or.
     &         IntTyp .eq. 4) Then
C The last argument is added to reformat the basis set to satisfy
C the requirement for the wave function analysis (see a2proc). The
C reformatting routine read a tmp file, GENB.TMP, created in the
C first MkVmol call and write the ZMAT.BAS file in the new format
C and the MOL file is always created from the the ZMAT.BAS (including
C the very first run). Tom Watson Jr. and Ajith Perera; 01/2011.
C
          Call MkVMOL (q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &         GenBy, Z(N11), IStat, BasNam)
C
         if (iflags(45).eq.2) then
            call mknddo(q, PGrp, NAtoms, NUnique, ZSym, IAtNum, GenBy,
     &                  Z(N11), IStat, BasNam)
         end if
      Else If (IntTyp .eq. 3) Then
cYAU         Call MkCadp (PGrp, NAtoms, NUnique, ZSym, IAtNum, Q,
cYAU     &                GenBy, IStat)
           write(*,*) '@GEOPT: MkCadp was temporarily removed.'
           write(*,*) '        This calculation cannot complete.'
           call errex
      End If
C
C************************ For the Future additions *****************
C
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
c THE GRADIENT IS IN Z(N7)
c THE HESSIAN  IS IN Z(N8)
C
      NX = MAX(TOTREDNCO, 3*NATOMS)
      IF (.NOT.OPTRES) CALL ARCHIVE(E,Z(N7),Z(N8),Z(N1),1,Z(N18))
      IF (IPRNT.EQ.999) WRITE(LUOUT,8888) (I,Z(I),I=1,N27)
      CALL PUTREC(20,'JOBARC','JODAOUT ',IONE,IJUNK)

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

      return
      END

