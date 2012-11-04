      SUBROUTINE MKabin(geom, PTGRP,NATMS,NUNIQUE,ZSYM,ATNR,GENBY,
     $     SCRATCH, ISTAT)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Purpose: Produce an VMOL .MOL file based on the ZMAT file
C
C Arguments:
C     PtGrp   Symbol for the point group symmetry to use in (input)
C             making the INP file.
C     NAtms   Number of atoms (input)
C     AtNr    Atomic numbers for each center (input)
C     NUnique Number of unique internal coordinates (input)
C     geom    Cartesian coordinates (input)
c             will not be the same as Q in FINDIF
C     GenBy   To hold the Genby list (internal)
C     IStat   Error codes (output)
C             = 0  Successful
C             = 1  Fatal error from dependent
C             = 3  I/O Error on AbiFil
C             = 5  I/O Error on ZFil
C             = 7  No non-dummy centers!
C             = 8  ZFil ends after Z-matrix
C             = 10 ZFil ends after JODA control info
C
C Common blocks used:
C     LUnits  Sets unit names & numbers
C     Coord
C
C Internal Variables:
C     BasNam  Holds the basis set name given in the ZMAT file
C     JobTtl  Holds the job title given in the ZMAT file
C     Wrk     Scratch string to copy lines form ZMAT to ABI
C
C Dependents
C     PitPtG  Translates a schoenflies symbol into Pitzer's form
C     WrPBas  Lookup & write basis set in Pitzer style
C
C Files used:
C     ZFil   Holds the basis specification & the "rest" of INP
C     AbiFil Gets written out here (kills previous one)
C     BasFil The basis set library
C
C Limitations:
C     Only looks up bases for the geometrically unique atoms.
C
C     Extra blank lines at the end of a "truncated" ZFil will cause
C     it to act like it encountered a real error on reading the file.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
      Character*(*) PtGrp, ZSym(NAtms)
      Integer AtNr(NAtms), GenBy(NAtms), IStat, AtoI,IZLoc(MxAtms)
      Character*80 JobTtl, Wrk
      Character*3 Ka(3)
      Character*6 Slask,NullSt
      Character*1 Lb
      Character*1 Indxx(150)

csb 1/97 Hold long basis set names
      integer baslen
      parameter (baslen=80)
      Character*(baslen) BasNam(MxAtms),BlnkBN,ScrBas,BasName

      Parameter (TLA = 1.D-9)
      INTEGER SHELLANG,SHELLLOC,SHELLSIZ,SHELLPRM
      INTEGER SHELLORB,SHOFFSET,PRIMORBT,PROFFSET
      integer fnblnk 
      double precision geom(3*mxatms)
C
      Logical Opn,bad123, GENBAS
      Integer OldLu
      Dimension Nord(2*MxAtms),Scratch(NAtms)
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
      COMMON /MOLSTR1/ SHELLANG(100),SHELLLOC(100),SHELLSIZ(100),
     &                 SHELLPRM(100),BASISEXP(10000),BASISCNT(10000)
      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IIATOM
      COMMON /MOLSTR3/ SHELLORB(100),SHOFFSET(100),PRIMORBT(100),
     &                 PROFFSET(100)
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /FLAGS/IFLAGS(100),IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /COORD/ Q(3*MXATMS), R(3*MXATMS), NCON(3*MXATMS),
     1     NR(MXATMS),ISQUASH(3*MXATMS),IATNUM(MXATMS),ATMASS(MXATMS),
     2     IUNIQUE(3*MXATMS),NEQ(3*MXATMS),IEQUIV(3*MXATMS,3*MXATMS),
     3     NOPTI(3*MXATMS), NATOMS
#include "io_units.par"
C
C     Length of BasNam is fixed by the basis library
C
C
C     This length is fixed by ABINITIO's input section
C
      Data Indxx /'1','2','3','4','5','6','7','8','9','0','A'
     &,'B','C','D','E','F','G','H','I','J','K','L','M','N'
     &,'O','P','Q','R','S','T','U','V','W','X','Y','Z','a'
     &,'b','c','d','e','f','g','h','i','j','k','l','m','n'
     &,'o','p','q','r','s','t','u','v','w','x','y','z',88*'*'/
csb      DATA BASIS /'PBS     ','STO-3G  ','DZ      ',
csb     &            'VDZ     ','3-21G   ','6-31G   ',
csb     &            '4-31G   ','TZ      ','PVTZ    ',
csb     &            '6-311G  ','DZP     ','PVDZ    ',
csb     &            'TZP     ','PVQZ    ','TZ2P    ',
csb     &            'WMR     ','6-31G*  ','6-31G** ',
csb     &            '6-311G* ','6-311G**','PV5Z    '/
C
C Initialize some things.
C
csb
      do i=1,baslen
        BlnkBN(i:i)=' '
      end do
      ScrBas=BlnkBN

      bad123 = .false.
      Slask='INTGRL'
      NullSt='      '
      LB='#'
      J=0
      I=0
      IGENER=1
      IF10=0
      IF18=1
      ID3=0
      IONE=1
      IRSTRT=0
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,INEWFD)
      IECP=IFLAGS(71)
      IECPU=0
      IPRINT=0
      INUTMP=0
      ISTATS=0
      ISTAT=0
      IOFFPRIM=1
      IOFFSHEL=1
C
      ICOUNT=0
C
C SEE IF A GENBAS FILE IS PRESENT
C
      INQUIRE(FILE='GENBAS',EXIST=GENBAS)
C
C     Cautiously open the ZMAT file
C
      Inquire (FILE=ZFil, OPENED=Opn)
      If (.not.Opn) Open (LuZ, FILE=ZFil, STATUS='OLD')
      Rewind LuZ
C
C     Read ZMAT: Get title, skip Z-matrix, skip parameters
C     NOTE: blank line after Z-matrix, and parameters required.
C
      Read (LuZ, '(A)', ERR=8500, END=8510) JobTtl
CCH   READ COORDINATES OR Z-MATRIX
      Do 200 i = 1, NAtms
         Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
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
      Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
      Do 205 i = 1, NUnique
         Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
 205  Continue
      Read (LuZ, '(A)', ERR=8500, END=8800) Wrk
C
C     Get the optimizer settings & another blank line after
C
      Call GtFlgs(0,IErr,IPrnt,INR,IContl,IVec,IDie,ICurvy,IStCrt,
     & IMxStp,IVib,IRecal,IntTyp,IDisFD,IGrdFD,ICNTYP,ISym,IBasis,
     & IDFGHI,BasName)
      If(IErr.eq.2)Goto 8800
      If(IErr.eq.1)Goto 8500
C      Read (LuZ, '(11I3)', ERR=8500, END=8800) IPrnt, INR, IContl,
C     $ IVec, IDie, ICurvy, IStCrt, IMxStp, IVib, IRecal, IntTyp
      If(IContl.eq.0)IContl = 4
      IF(IPrnt.ge.100)then
      Write (LuOut, '(4X,A,I3)') 'Print Level: ', IPrnt
      Write (LuOut, '(4X,A,I3)') 'Step-taking: ', INR
      Write (LuOut, '(4X,A,I3)') 'E-Vec to follow: ',IVec
      Write (LuOut, '(4X,A,I3)') 'Die: ',IDie
      Write (LuOut, '(4X,A,I3)') 'Use curvilinear xform: ',ICurvy
      Write (LuOut, '(4X,A,I3)') 'Max Step: ',IMxStp
      Write (LuOut, '(4X,A,I3)') 'Step Contol: ', IStCrt
      Write (LuOut, '(4X,A,I3)') 'Vibrational Analysis: ', IVib
      Write (LuOut, '(4X,A,I3)') 'Recalculation of Hessian: ',
     &IRecal
      Write (LuOut, '(4X,A,I3)') 'Convergence cutoff: ',IContl
      Write (LuOut, '(4X,A,I3)') 'Integral type: ',IntTyp
      Endif
      Read (LuZ, '(A)', ERR=8500, END=8810) Wrk
C
C     Now read the basis sets - each line corresponds to a (non-dummy)
C     center in the Z-matrix, in the order of appearance.  We put this
C     into an array length NAtms, with blanks for dummies.  This makes
C     it easier to match the basis with the right coordinates later.
C     NOTE: blank line after basis names.
C
      Do 250 i = 1, NAtms
         BasNam(i) = BlnkBN
         If (AtNr(i) .ne. 0 .and. IBasis .eq. 0) then
            Read (LuZ, '(A)', ERR=8500, END=8510, IOSTAT=IOS) BasNam(i)
         ElseIf(AtNr(i) .ne. 0 .and. IBasis .ne. 0)then
            BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
     &          BasName(1:linblnk(BasName))
csb            izz=index(Basis(IBasis),'*')
csb            ixx=index(Basis(IBasis),'**')
csb            ipop=0
csb            if(ixx.ne.0.and.atmass(i).gt.4.5)ipop=1
csb            if(izz.ne.ixx.and.atmass(i).lt.4.5)ipop=1
csb            BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
csb     &                Basis(IBasis)(1:linblnk(Basis(IBasis))-ipop)
         EndIf
 250  Continue
      If(IBasis.eq.0)Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
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
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Let's find out which centers are equivalent
C
      Call SymEqv ( NAtms, GenBy)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Count number of centers having *different* atomic numbers in this grou
C   Top half of Nord contains the number of symmetry equivalent atoms
C   for a given atomic number.
C
      Do 312 I=1,NAtms
       If(GenBy(i).eq.0.and.AtMass(i).gt.1.D-2)Then
        J=J+1
        Scratch(J)=AtMass(I)
        Nord(J)=I
       EndIf
 312  Continue
      If(IVIB.ne.2)then
c       Call PikSr2(J,Scratch,Nord)
       NonTyp=1
       If(J.Eq.1)Then
        Nord(MxAtms+1)=1
       Else
        Do 313 I=2,J
         Dif=Dabs(Scratch(i)-Scratch(i-1))
         If(NonTyp.eq.1)IBot=0
c         If(Dif.lt.1.D-14)GoTo 314
         Nord(MxAtms+Nontyp)=I-1-IBot
         IBot=Nord(MxAtms+Nontyp)+IBot
         NonTyp=NonTyp+1
  314    If(I.eq.J)Nord(MxAtms+Nontyp)=I-Ibot
 313    Continue
       EndIf
C
C For finite difference calcs, reset this so that it writes out the
C  atoms to the MOL
C  file in exactly the same order as they appear in the ZMAT.  This
C  ordering is assumed in FINDIF!
C
      ElseIf(ivib.eq.2)then
       NonTyp=J
       Do 315 I=1,J
        Nord(MxAtms+I)=1
315    Continue
      EndIf
C
C Now get symmetry information written by SYMMETRY to file VMINF.
C
      Opn = .False.
      OldLu = 0
      IF(PTGRP.EQ.'C1 ')THEN
       INQUIRE(FILE='VMLSYM',EXIST=Opn)
       IF(OPN)THEN
        OPEN(UNIT=30,FILE='VMLSYM',FORM='UNFORMATTED',
     &       STATUS='OLD')
        CLOSE(UNIT=30,STATUS='DELETE')
       ENDIF
       GOTO 5400
      ENDIF
      Inquire (FILE='VMLSYM', OPENED=Opn, NUMBER=OldLu)
      If (Opn) Close (OldLu)
C
      Open (30, FILE='VMLSYM', STATUS='OLD',FORM='UNFORMATTED')
      Rewind(30)
      Read (30,Err=9400)NSymOp,(KA(i),i=1,3)
      Close(30, Status='Delete')
5400  If(PtGrp.EQ.'C1 ')Then
       Do 412 I=1,3
        KA(I)='   '
412    Continue
       NSymOp=0
      EndIf
C
C     Now we've got (pretty much) all the information we need, so
C     we can start writing out the MOL file.
C     Begin with a very cautious opening of the file.
C
      Opn = .False.
      Inquire (file=abifil, OPENED=Opn)
      If ( .not. Opn ) Open (Luabi, FILE=abiFil, FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      Rewind Luabi
C     

       IF(GENBAS)THEN
          CALL MKARGOS(KA,NSYMOP,PTGRP,NONTYP,ZSYM,BASNAM,NORD,GEOM,
     $           JOBTTL,IATNUM,luabi,natms,izloc)
       else
          write(luerr,9913)
 9913     format(' @MKDECK-F, Sorry, GENBAS file required')
          call errex
       endif
C
C     Dump vector to disk which relates the atoms listed in the
C     MOL file to their position in the user supplied Z-matrix.
C
      CALL PUTREC(20,'JOBARC','MAP2ZMAT',NATMS,IZLOC)
C
C This better not happen.
C
      If (Nontyp .eq. 0) then
         IStat = 7
         Write (LuErr, 9970)
         Close (LuZ)
         Close (LuAbi)
         Return
 9970    Format (' @MKDECK-F, No non-dummy centers in Z-matrix.')
      EndIf
C
C     The rest of this just copied from ZMAT
C
C
C     Make sure every body gets closed & stuff
C
 8000 Close (LuAbi)
      Close (LuZ)
C
C     WRITE SCF INPUT FILE
C
      IF(IFLAGS(17).NE.999) CALL SCFINP(JOBTTL,LUABI)
C
C WRITE BASIS SET INFORMATION TO JOBARC
C
      IONE=1
      NTOTPRIM=IOFFPRIM-1
      NTOTSHEL=IOFFSHEL-1
      CALL PUTREC(20,'JOBARC','ATOMCHRG',NATOMS,IATNUM)
      IF(INEWFD.EQ.0)THEN
       CALL PUTREC(20,'JOBARC','NTOTSHEL',IONE,NTOTSHEL)
       CALL PUTREC(20,'JOBARC','NTOTPRIM',IONE,NTOTPRIM)
       CALL PUTREC(20,'JOBARC','BASISEXP',NTOTPRIM*IINTFP,BASISEXP)
       CALL PUTREC(20,'JOBARC','BASISCNT',NTOTPRIM*IINTFP,BASISCNT)
       CALL PUTREC(20,'JOBARC','SHELLSIZ',NTOTSHEL,SHELLSIZ)
       CALL PUTREC(20,'JOBARC','SHELLPRM',NTOTSHEL,SHELLPRM)
       CALL PUTREC(20,'JOBARC','SHELLANG',NTOTSHEL,SHELLANG)
       CALL PUTREC(20,'JOBARC','SHELLLOC',NTOTSHEL,SHELLLOC)
       CALL PUTREC(20,'JOBARC','SHOFFSET',NONTYP,SHOFFSET)
       CALL PUTREC(20,'JOBARC','SHELLORB',NONTYP,SHELLORB)
       CALL PUTREC(20,'JOBARC','PROFFSET',NONTYP,PROFFSET)
       CALL PUTREC(20,'JOBARC','PRIMORBT',NONTYP,PRIMORBT)
      ENDIF
      Return
C
C     Get here via an I/O error on VMLSYM
C
 9400 IStat=3
      Write (LuErr, 9930)'VMLSYM'
      Close(LuAbi)
      Close(LuZ)
      Return
C
C     Get here via an I/O error on LuZ
C
 8500 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9950) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9950) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 8510 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9955) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9955) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
C
C     If the ZFil ends after the Z-matrix or the optimization control
C     info, it probably means that JODA is being used "stand-alone"
C     and we want to allow for that possibility.
C
 8800 IStat = 8
      If(IntTyp.ne.0)then
       Write (LuOut, 9800) ZFil, MolFil
       Close (LuVMol)
      else
       Write (LuOut, 9800) ZFil, AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 8810 IStat = 10
      If(IntTyp.ne.0)then
       Write (LuOut, 9810) ZFil, MolFil
       Close (LuVMol)
      else
       Write (LuOut, 9810) ZFil, AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 9930 Format (' @MKDECK-F, I/O error on file ',A,'.')
 9950 Format (' @MKDECK-F, I/O error on file ',A,'.')
 9955 Format (' @MKDECK-F, Premature end-of-file on ',A,'.')
 9800 Format (/1X,A,' ends after Z-matrix - cannot ',
     $   'finish making ',A,'.')
 9810 Format (/1X,A,' ends after JODA control info - cannot ',
     $   'finish making ',A,'.')
      End
