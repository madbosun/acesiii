      Subroutine MkCadp (PtGrp, NAtms, NUnique, ZSym, AtNr, QX,
     &GenBy, IStat)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Purpose: Produce a CADPAC input file based on the ZMAT file
C
C Arguments:
C     PtGrp   Symbol for the point group symmetry to use in (input)
C             making the INP file.
C     NAtms   Number of atoms (input)
C     AtNr    Atomic numbers for each center (input)
C     NUnique Number of unique internal coordinates (input)
C     Q       Cartesian coordinates (input)
C     GenBy   To hold the Genby list (internal)
C     IStat   Error codes (output)
C             = 0  Successful
C             = 1  Fatal error from dependent
C             = 3  I/O Error on CadFil
C             = 5  I/O Error on ZFil
C             = 7  No non-dummy centers!
C             = 8  ZFil ends after Z-matrix
C             = 10 ZFil ends after JODA control info
C
C Common blocks used:
C     LUnits  Sets unit names & numbers
C
C Internal Variables:
C     BasNam  Holds the basis set name given in the ZMAT file
C     JobTtl  Holds the job title given in the ZMAT file
C     Wrk     Scratch string to copy lines form ZMAT to ABI
C     PitSym  Holds the point group label in Pitzer's form
C     PitOrd  Order of the axis - for Pitzer
C     IOS     Gets the status code from all I/O statements
C     ISect   Counts sections in ZMAT and ABI files to help
C             locate errors.
C
C Dependents
C     PitPtG  Translates a schoenflies symbol into Pitzer's form
C     WrPBas  Lookup & write basis set in Pitzer style
C
C Files used:
C     ZFil   Holds the basis specification & the "rest" of INP
C     CadFil Gets written out here (kills previous one)
C     BasFil The basis set library
C
C Limitations:
C     Only looks up bases for the geometrically unique atoms.
C
C     Extra blank lines at the end of a "truncated" ZFil will cause
C     it to act like it encountered a real error on reading the file.
C
C Diagnostics
C     Printed for most I/O errors and end-of-file conditions.  The
C     IOSTAT result code is reported where appropriate, along with
C     the section of the file in which the error occurred.  The
C     IOSTAT code is system dependent and probably requires manuals
C     for the machine.
C
C     ZFIL file sections:
C     1) Title and Z-matrix
C     2) Optimization control input
C     3) Basis sets
C     4) Point group input
C     5) Remainder of file
C
C     CadFil file sections:
C     1) Title
C     2) Point group
C     3) Atom label & coordinates
C     4) Basis
C     5) Remainder of file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
      Character*(*) PtGrp, ZSym(NAtms)
      Integer AtNr(NAtms), GenBy(NAtms), IStat, IOS, ISect
      Integer AtoI
      Character*80 JobTtl, Wrk
      Character*5 PitSym
      Integer PitOrd

csb 1/97 Hold long basis set names
      integer baslen
      parameter (baslen=80)
      character*(baslen) BasNam(MxAtms),BlnkBN,ScrBas,BasName

      Logical Opn,bad123
      Integer OldLu
      Dimension QX(3*NAtms)
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
 
#include "io_units.par"
 
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      Common /SymInf/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      COMMON /COORD/ Q(3*MXATMS), R(3*MXATMS), NCON(3*MXATMS),
     1     NR(MXATMS),ISQUASH(3*MXATMS),IATNUM(MXATMS),ATMASS(MXATMS),
     2     IUNIQUE(3*MXATMS),NEQ(3*MXATMS),IEQUIV(3*MXATMS,3*MXATMS),
     3     NOPTI(3*MXATMS), NATOMS
C
C     Length of BasNam is fixed by the basis library
C
C
C     This length is fixed by ABINITIO's input section
C
c$$$      DATA BASIS /'PBS     ','STO-3G  ','DZ      ',
c$$$     &            'VDZ     ','3-21G   ','6-31G   ',
c$$$     &            '4-31G   ','TZ      ','VTZ     ',
c$$$     &            '6-311G  ','DZP     ','VDZP    ',
c$$$     &            'TZP     ','VTZP    ','TZ2P    ',
c$$$     &            'WMR     ','6-31G*  ','6-31G** ',
c$$$     &            '6-311G* ','6-311G**','PV5Z    '/
C
C     Cautiously open the ZMAT file
C
      PitOrd = 0
      Opn = .False.
      OldLu = 0
      Inquire (FILE=ZFil, OPENED=Opn, NUMBER=OldLu)
      If (Opn) Close (OldLu)
C
      Open (LuZ, FILE=ZFil, STATUS='OLD')
      Rewind LuZ
C
C     Read ZMAT: Get title, skip Z-matrix, skip parameters
C     NOTE: blank line after Z-matrix, and parameters required.
C
      bad123 = .false.
      ISect = 1
      Read (LuZ, '(A)', ERR=8500, END=8510, IOSTAT=IOS) JobTtl
      Do 200 i = 1, NAtms
         Read (LuZ,'(A)', ERR=8500, END=8510, IOSTAT = IOS) Wrk
C
C     ADDITION (JFS, 4/1990): Check third atom in Z-matrix in
C     order to detect 2--1--3 connectivity.
C
         If(I.Eq.3)then
          iblank=index(wrk,' ')
          istrt2=fnblnk(wrk(iblank+1:))+iblank
          iendd2=index(wrk(istrt2+1:),' ')+istrt2-1
          izz=AtoI(wrk(istrt2:iendd2))
          if(izz.eq.1.and.zsym(1)(1:1).ne.'X'.and.zsym(2)(1:1).ne.'X')
     &       bad123=.true.
         endif
C     END OF ADDITION
 200  Continue
      Read (LuZ,'(A)', ERR=8500, END=8510, IOSTAT=IOS) Wrk
      Do 205 i = 1, NUnique
         Read (LuZ, '(A)', ERR=8500, END=8510, IOSTAT=IOS) Wrk
 205  Continue
      Read (LuZ, '(A)', ERR=8500, END=8800, IOSTAT=IOS) Wrk
C
C     Get the optimizer settings & another blank line after
C
      ISect = 2
      Call GtFlgs(0,IERR,IPRNT,INR,ICONTL,IVEC,IDIE,ICURVY,ISTCRT,
     &IMXSTP,IVIB,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     &IDFGHI,BasName)
      IF(IERR.EQ.1)GOTO 8500
      IF(IERR.EQ.2)GOTO 8800
C      Read (LuZ, '(10I3)', ERR=8500, END=8800) IPrnt, INR, IContl,
C     $ IVec, IDie, ICurvy, IStCrt, IMxStp, IVib, IRecal
      If(IContl.eq.0)IContl = 4
      IF(IPrnt.ge.100)then
      Write (LuOut, '(1X,A,I3)') 'Print Level: ', IPrnt
      Write (LuOut, '(1X,A,I3)') 'Step-taking: ', INR
      Write (LuOut, '(1X,A,I3)') 'E-Vec to follow: ',IVec
      Write (LuOut, '(1X,A,I3)') 'Die: ',IDie
      Write (LuOut, '(1X,A,I3)') 'Use curvilinear xform: ',ICurvy
      Write (LuOut, '(1X,A,I3)') 'Max Step: ',IMxStp
      Write (LuOut, '(1X,A,I3)') 'Step Contol: ', IStCrt
      Write (LuOut, '(1X,A,I3)') 'Vibrational Analysis: ', IVib
      Write (LuOut, '(1X,A,I3)') 'Recalculation of Hessian: ',
     &IRecal
      Write (LuOut, '(1X,A,I3)') 'Convergence cutoff: ',IContl
      Endif
      Read (LuZ, '(A)', ERR=8500, END=8810, IOSTAT=IOS) Wrk
C
C     Now read the basis sets - each line corresponds to a (non-dummy)
C     center in the Z-matrix, in the order of appearance.  We put this
C     into an array length NAtms, with blanks for dummies.  This makes
C     it easier to match the basis with the right coordinates later.
C     NOTE: blank line after basis names.
C
      ISect = 3
      Do 250 i = 1, NAtms
         BasNam(i) = BlnkBN
         If (AtNr(i) .ne. 0 .and. IBasis .eq. 0) then
            Read (LuZ, '(A)', ERR=8500, END=8510, IOSTAT=IOS) BasNam(i)
         ElseIf(AtNr(i) .ne. 0 .and. IBasis .ne. 0)then
            BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
     &          BasName(1:linblnk(BasName))
c$$$            izz=index(Basis(IBasis),'*')
c$$$            ixx=index(Basis(IBasis),'**')
c$$$            ipop=0
c$$$            if(ixx.ne.0.and.atmass(i).gt.4.5)ipop=1
c$$$            if(izz.ne.ixx.and.atmass(i).lt.4.5)ipop=1
c$$$            BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
c$$$     &                Basis(IBasis)(1:linblnk(Basis(IBasis))-ipop)
         EndIf
 250  Continue
      If(IBasis.eq.0)
     &Read (LuZ, '(A)', ERR=8500, END=8510, IOSTAT=IOS) Wrk
C ADDITION (JFS,4/90)
C Now take care of situation which occurs when 2--1--3 Z-matrix
C  specification and nonstandard basis set input is used.  In this case,
C  the basis sets for atoms #1 and #2 must be switched.
C
      If(Bad123.and.IBasis.eq.0)then
       ScrBas=BasNam(1)
       BasNam(1)=BasNam(2)
       BasNam(2)=ScrBas
      EndIf
C END OF ADDITION
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Let's find out which centers are equivalent.  This call is
C     skipped in finite difference calculations.
C
      If(IVib.ne.2)
     &Call SymEqv ( NAtms, GenBy)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Now we've got (pretty much) all the information we need, so
C     we can start writing out the INP file.
C     Begin with a very cautious opening of the file.
C
      Inquire (FILE=CadFil, OPENED=Opn, NUMBER=OldLu)
      If (Opn) Close (OldLu)
C
      Open (LuCad, FILE=CadFil, FORM='FORMATTED', STATUS='UNKNOWN')
      Rewind LuCad
C
C     This we got from the beginning of the ZMAT file
C
      ISect = 1
      Write (LuCad, '(A)', ERR=8300, IOSTAT=IOS) 'TITLE'
      Write (LuCad, '(A)', ERR=8300, IOSTAT=IOS)
     $   JobTtl(:linblnk(JobTtl))
C
C     Now for the point group symbol
C
      Call PitPtG (PtGrp, PitSym, PitOrd)
      ISect = 2
      Write (LuCad, 9010, ERR=8300, IOSTAT=IOs) PitSym, PitOrd
 9010 Format ('SYMMETRY'/A5,I5/'END')
C
C     Now to write out the geometry/basis information
C
      Write (LuCad, *, ERR=8300, IOSTAT=IOS)
      Write (LuCad, '(A)', ERR=8300, IOSTAT=IOS) 'ATOMS'
      IPit = 0
      Do 100 i = 1, NAtms
         If (AtNr(i) .ne. 0 .AND. GenBy(i) .eq. 0) then
            IPit = IPit + 1
            ISect = 3
            Write (LuCad, 9000, ERR=8300, IOSTAT=IOS) ZSym(i),
     $         Float(AtNr(i)), (QX(3*(i-1)+j),j=1,3)
 9000       Format (A5, 3X, 2X, F5.0, 3F20.12)
C
            ISect = 4
            Call WrPBas (LuBasL,BasFil,BasNam(i),LuCad,1,2,IStat)
            Write (LuCad, '(A)', ERR=8300, IOSTAT=IOS) 'END'
C
C           Handle errors in WrPBas
C
            If (Mod(IStat,2) .eq. 1) then
               IStat = 1
               Write (LuErr, 9910)
               Close (LuZ)
               Close (LuCad)
               Return
            ElseIf (IStat .eq. 8) then
               IStat = 3
               Write (LuErr, 9935) CadFil, LuCad
               Close (LuZ)
               Close (LuCad)
               Return
            EndIf
         EndIf
 100  Continue
 9910 Format (' @MKCADP-F, Dependent terminated with FATAL error.')
 9935 Format (' @MKCADP-F, Dependent reports error accessing file ',
     $   A,' on unit ',I3,'.')
C
C     This really better not happen...
C
      If (Ipit .eq. 0) then
         IStat = 7
         Write (LuErr, 9970)
         Close (LuZ)
         Close (LuCad)
         Return
 9970    Format (' @MKCADP-F, No non-dummy centers in Z-matrix.')
      EndIf
C
C     Finish off the "atoms" section
C
      Write (LuCad, '(A)', ERR=8300, IOSTAT=IOS) 'END'
C
C     The rest of this just copied from ZMAT
C
      Write (LuCad, * , ERR=8300, IOSTAT=IOS)
 300  Continue
      ISect = 5
      Read (LuZ, '(A)', END=500, IOSTAT=IOS) Wrk
      Write (LuCad, '(A)', ERR=8300, IOSTAT=IOS) Wrk(:linblnk(Wrk))
      Goto 300
 500  Continue
C
C     Make sure every body gets closed & stuff
C
 8000 Close (LuCad)
      Close (LuZ)
      Return
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Get here via an I/O error on LuCad
C
 8300 IStat = 3
      Write (LuErr, 9930) IOS, CadFil, ISect
      Close (LuCad)
      Close (LuZ)
      Return
 9930 Format (' @MKCADP-F, I/O error ',I5,' on file ',A,' in section ',
     $   I2,'.')
C
C     Get here via an I/O error on LuZ
C
 8500 IStat = 5
      Write (LuErr, 9950) IOS, ZFil, ISect
      Close (LuCad)
      Close (LuZ)
      Return
 8510 IStat = 5
      Write (LuErr, 9955) ZFil, ISect
      Close (LuCad)
      Close (LuZ)
      Return
C
C     If the ZFil ends after the Z-matrix or the optimization control
C     info, it probably means that JODA is being used "stand-alone"
C     and we want to allow for that possibility.
C
 8800 IStat = 8
      Write (LuOut, 9800) ZFil, CadFil
      Close (LuCad)
      Close (LuZ)
      Return
 8810 IStat = 10
      Write (LuOut, 9810) ZFil, CadFil
      Close (LuCad)
      Close (LuZ)
      Return
 9950 Format (' @MKCADP-F, I/O error ',I5,' on file ',A,' in section ',
     $   I2,'.')
 9955 Format (' @MKCADP-F, Premature end-of-file on ',A,' in section ',
     $   I2,'.')
 9800 Format (/1X,A,' ends after Z-matrix - cannot ',
     $   'finish making ',A,'.')
 9810 Format (/1X,A,' ends after JODA control info - cannot ',
     $   'finish making ',A,'.')
      End
