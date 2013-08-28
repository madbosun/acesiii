      SUBROUTINE FSPUTT1(Z,ILIST1,ILIST2,IOTYPE,SYMTYP)
C
C THIS ROUTINE READS A "STANDARD" TWO INDEX LIST IN A
C NONSTANDARD WAY.  IT ALLOWS ONE TO READ LISTS OF THE TYPE:
C
C           W(PQ)  WITH IOTYPE = 'N1'
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
      DOUBLE PRECISION Z(1)
      CHARACTER*2 IOTYPE
      COMMON // ICORE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FSIOLOC/ IFSPOS(8,22,8),ISCRLC,IFSPOSZ(8,22,8)
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
      IOFF=1
C
C     Note that T1 lists are always Irrep=1 for symmetry purposes, 
C     though they may not be stored on the first sublist.
C
      If (IOType(1:2) .eq. 'NN') then
         LType = 0
      Else
         LType = iszeq(8, FSDsLbl, 1, IOType(1:2))
         If (LType .ne. 0) then
            IPosLeft = IFSPos(  1, SymTyp, LType)
            IZPosL   = IFSPosZ( 1, SymTyp, LType)
            NSize    = FSDPD(   1, SymTyp, LType)
         Else
            LType = iszeq(TMxDsTy, TDsLbl, 1, IOType(1:2))
            If (LType .ne. 0) then
               IPosLeft = TIOPos (1, SymTyp, LType)
               IZPosL   = TIOZPos(1, SymTyp, LType)
               NSize    = TIODPD (1, SymTyp, LType)
            Else
               Write (6, 9000) IOType(1:2)
               Call ErrEx
            EndIf
         EndIf
      EndIf
C
 9000 Format(T3,'@FSPUTT1-F, Unknown index type ', A, ' requested.')
C
C     Allow 'NN' requests to work too
C
      If (LType .eq. 0) then
         CALL PUTLST(Z,1,1,1,ILIST1,ILIST2)
      Else
         If (IZPosL .ne. 0) then
            Do 15 I = 1, NSize
               Z(I) = Z(I) * ICore(IZPosL + I - 1) 
 15         Continue
         EndIf
         IF(TIOFULL.EQ.1) CALL GETLST(ICORE(ISCRLC),1,1,1,ILIST1,ILIST2)
         CALL SCATTER(NSIZE,ICORE(ISCRLC),ICORE(IPOSLEFT),Z)
         CALL PUTLST(ICORE(ISCRLC),1,1,1,ILIST1,ILIST2)
      EndIf
C
      RETURN
      END
