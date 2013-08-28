      SUBROUTINE FSPUT(Z,ISTART,NUMPUT,IDUMMY,IRREP,LIST,IOTYPE)
C
C THIS ROUTINE WRITES A "STANDARD" INTEGRAL OR AMPLITUDE LIST IN A
C NONSTANDARD WAY.  IT ALLOWS ONE TO READ LISTS OF THE TYPE:
C
C           W(PQ,R1)  WITH IOTYPE = 'NNNA'
C
C WHERE 1 IS A SPECIAL TYPE OF ORBITAL (ACTIVE OR INACTIVE).  ANY
C INDEX CAN REFER TO A "SPECIAL" ORBITAL
C
C  'N' - NORMAL CASE.  READ ALL ORBITALS OF THIS TYPE
C  'A' - ACTIVE ORBITALS ONLY.
C  'I' - INACTIVE ORBITALS ONLY.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION Z(*)
      CHARACTER*4 IOTYPE
      COMMON // ICORE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FSIOLOC/ IFSPOS(8,22,8),ISCRLC, IFSPOSZ(8,22,8)
      COMMON /SYMPOP  / IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
C
C     It is convenient to consider all of the FSDPD__ arrays as a
C     single large array.
C
      Integer FSDPD(8,22,8)
      Equivalence (FSDPDAN(1,1), FSDPD(1, 1, 1) )
C
C     Common block to hold labels for the index types
C
      Character*2 FSDsLbl(8)
      Common /FSIOC/ FSDsLbl
C
C&      Include 'tempio.inc'
C
C     For temporary I/O vectors, etc.
C
C     TIOBas  Pointer to space allocated in ICore for temporary I/O
C     TIOSiz  Size of space allocated in ICore for temporary I/O
C     TIOPos  Pointer to gather vector in ICore
C     TIOZPos Position of zero list in ICore
C     TIODPD  Distribution size
C     TIOFull If it is set to 1, the lists are 'right' after 
C             FSPUT and FSPUTT1. Otherwise only the active part
C	      can be used.
C     TIOBufSz FSIO buffer size -- needed by FSMAP to properly
C              handle writing lists of TIOFull=1.
C
C     Dimensions are: 8 irrpes, 22 list types (even though we
C     don't use most of them), up to 5 different types of distributions.
C
C     TIOZPos is only used by UPLE and UPLT packing schemes.  For the
C     others, TIOZPos(*, IType, IDsTyp) must be 0 or FSGET/FSPUT will be
C     very confused.
C
      Integer TMxDsTy
      Parameter (TMxDsTy = 5)
C
      Integer TIOBas, TIOSiz, TIOPos(8, 22, TMxDsTy),
     $   TIOZPos(8, 22, TMxDsTy), TIODPD(8, 22, TMxDsTy), TIOFull,
     $   TIOBufSz
C
C     Text labels for the names of the temp. I/O distributions
C     (slowest dimension of TIO info arrays).
C
      Character*2 TDsLbl(TMxDsTy)
C
      Common /TempIO/ TIOBas, TIOSiz, TIOPos, TIOZPos, TIODPD, TIOFull,
     $   TIOBufSz
      Common /TempIOC/ TDsLbl
C
C LOCAL VARIABLES
C
C     LType   Index type of left-hand index pair.  0 is reserved for
C             'NN', which means no gather is done; otherwise, a positive
C             integer corresponding to one of the standard index types
C             or one of the temporary IO types (no distinction is
C             required once IPosLeft, NSize, and IZPosL are set).
C     RType   Same as LType for the right-hand index pair.
C     IPosLeft Pointer to gather vector for LHS indices.
C     IPosRght Pointer to gather vector for RHS indices.
C     IZPosL   Pointer to "zero list" vector for LHS indices. 0 if
C              no extra processing is required.
C     IZPosR   Same as IZPosL for RHS.
C     I        Loop index
C     IFactor  A particular factor from the "zero list" for the ket
c              indices.
C
      Integer LType, RType, IPosLeft, IPosRght, IZPosL, IZPosR, I,
     $   IFactor
C
      IOFF=1
      IPosLeft = 0
      IPosRght = 0
      IZPosL = 0
      IZPosR = 0
C
      SYTYPL=ISYTYP(1,LIST)
      SYTYPR=ISYTYP(2,LIST)
C
C DISTINGUISH LEFT AND RIGHT TYPES.  SET UP POINTER ADDRESSES.
C
      If (IOType(1:2) .eq. 'NN') then
         LType = 0
      Else
         LType = iszeq(8, FSDsLbl, 1, IOType(1:2))
         If (LType .ne. 0) then
            IPosLeft = IFSPos(  Irrep, SyTypL, LType)
            IZPosL   = IFSPosZ( Irrep, SyTypL, LType)
            NSize    = FSDPD(   Irrep, SyTypL, LType)
         Else
            LType = iszeq(TMxDsTy, TDsLbl, 1, IOType(1:2))
            If (LType .ne. 0) then
               IPosLeft = TIOPos (Irrep, SyTypL, LType)
               IZPosL   = TIOZPos(Irrep, SyTypL, LType)
               NSize    = TIODPD (Irrep, SyTypL, LType)
cSSS               Write (6, 8000) Irrep, List, 'Left', LType, SyTypL,
cSSS     $            IPosLeft, IZPosL
            Else
               Write (6, 9000) IOType(1:2), 'left'
               Call ErrEx
            EndIf
         EndIf
      EndIf
C
      If (IOType(3:4) .eq. 'NN') then
         RType = 0
      Else
         RType = iszeq(8, FSDsLbl, 1, IOType(3:4))
         If (RType .ne. 0) then
            IPosRght = IFSPos(  Irrep, SyTypR, RType)
            IZPosR   = IFSPosZ( Irrep, SyTypR, RType)
         Else
            RType = iszeq(TMxDsTy, TDsLbl, 1, IOType(3:4))
            If (RType .ne. 0) then
               IPosRght = TIOPos (Irrep, SyTypR, RType)
               IZPosR   = TIOZPos(Irrep, SyTypR, RType)
cSSS               Write (6, 8000) Irrep, List, 'Right', RType, SyTypR,
cSSS     $            IPosRght, IZPosR
            Else
               Write (6, 9000) IOType(3:4), 'right'
               Call ErrEx
            EndIf
         EndIf
      EndIf
C
 8000 Format(1X, '@FSPUT-I, List [', I3, ',', I3, '] ', A5,
     $   ' indices are dist., ', I2, ' sym. type ', I2,/11X,
     $   'Gather vector at ', I6, ' zero vector at ', I6, '.')
 9000 Format(T3,'@FSPUT-F, Unknown index type ', A, ' requested ',
     $   'for ', A, ' index pair.')
C
C     Sanity check for gather vector pointers
C
      If (LType .ne. 0 .AND. IPosLeft .eq. 0) then
         Write (6, 9100) Irrep, List, IOType, 'left'
         Call ErrEx
      EndIf
      If (RType .ne. 0 .AND. IPosRght .eq. 0) then
         Write (6, 9100) Irrep, List, IOType, 'right'
         Call ErrEx
      EndIf
 9100 Format(1X, '@FSPUT-F, Bad gather vector pointer for list [', I3,
     $   ',', I3, '] ', A, 1X, A, ' indices.')
C
C CASE I.  ONLY BRA INDICES ARE NONSTANDARD
C
      IF( RType .eq. 0 )THEN
       DO 10 IPUT=ISTART,NUMPUT
C
C       Apply factors from "zero list" if necessary
C
        If (IZPosL .ne. 0) then
           Do 15 I = 0, NSize - 1
              Z(IOff + I) = ICore(IZPosL + I ) * Z(IOff + I)
 15        Continue
        EndIf
C
        IF(TIOFULL.EQ.1) CALL GETLST(ICORE(ISCRLC),IPUT,1,1,IRREP,LIST)
        CALL SCATTER(NSIZE,ICORE(ISCRLC),ICORE(IPOSLEFT),Z(IOFF))
        CALL PUTLST(ICORE(ISCRLC),IPUT,1,1,IRREP,LIST)
C
        IOFF=IOFF+NSIZE
10     CONTINUE
C
C CASE II. ONLY KET INDICES ARE NONSTANDARD
C
      ELSEIF( LType .eq. 0 )THEN
       NSIZE=IRPDPD(IRREP,SYTYPL)
       DO 20 IPUT=ISTART,NUMPUT
        IREALPUT=ICORE(IPOSRGHT+IPUT-1)
C
C       Do we have a factor from the "zero list" to account for?
C
        If (IZPosR .ne. 0) then
           IFactor = ICore(IZPosR + IPut - 1)
        Else
           IFactor = 1
        EndIf
C
C       Might have to scale everything
C
        If (IFactor .ne. 1 .AND. IFactor .ne. 0) then
           Do 25 I = 0, NSize - 1
              Z(IOff + I) =  Z(IOff + I) * IFactor 
 25        Continue
        EndIf
C
C       If factor is 0, don't even try to write it
C
        If (IFactor .ne. 0) then
           CALL PUTLST(Z(IOFF),IREALPUT,1,1,IRREP,LIST)
        EndIf
C
        IOFF=IOFF+NSIZE
20     CONTINUE
C
C CASE III. BOTH NONSTANDARD
C
      ELSE
       DO 30 IPUT=ISTART,NUMPUT
        IREALPUT=ICORE(IPOSRGHT+IPUT-1)
        If (IZPosR .ne. 0) then
           IFactor = ICore(IZPosR + IPut - 1)
        Else
           IFactor = 1
        EndIf
C
C       If factor is 0, don't even try to read/write it
C
        If (IFactor .ne. 0) then
C
C       Apply factors from "zero list" if necessary
C
        If (IZPosL .ne. 0) then
           Do 35 I = 0, NSize - 1
              Z(IOff + I) = Z(IOff + I) * (IFactor * ICore(IZPosL + I ))
 35        Continue
        EndIf
           IF(TIOFULL.EQ.1) 
     $        CALL GETLST(ICORE(ISCRLC),IREALPUT,1,1,IRREP,LIST)
           CALL SCATTER(NSIZE,ICORE(ISCRLC),ICORE(IPOSLEFT),Z(IOFF))
           CALL PUTLST(ICORE(ISCRLC),IREALPUT,1,1,IRREP,LIST)
        EndIf
C
        IOFF=IOFF+NSIZE
30     CONTINUE
C
      ENDIF
C
      RETURN
      END
