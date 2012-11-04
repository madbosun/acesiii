      SUBROUTINE FSIODRV(ICORE,MAXCOR,IUHF,IDO,IACTIVE,
     &                   NACTIVE,NBAS,I0,TFULL)
C
C DRIVER FOR CONSTRUCTION OF THE GATHER-SCATTER I/O VECTORS 
C NEEDED FOR CALCULATIONS INVOLVING ACTIVE OR INACTIVE INDICES
C
      IMPLICIT INTEGER (A-Z)
#include "maxbasfn.par"
      DATA MXBAS /MAXBASFN/
C
C     NOTE: ICore is really MaxCor long, but we can't use anything
C     below I0.
C
      DIMENSION ICORE(I0+MAXCOR),IDO(22),SPIN(2,22),NACTIVE(2)
      CHARACTER*3 TYPE(2,22),LFTTYP,RHTTYP
      CHARACTER*4 PACK(22),TYPEX
      DIMENSION IACTIVE(MAXBASFN,2),IVECTOR(MAXBASFN,3,2),
     &          POPLFT(8),POPRHT(8)
C
      Common /FSORBS/ IVector
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM   / POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      Integer NOccO(2), NVrtO(2)
      Common /Info/ NOccO, NVrtO
C
C COMMON BLOCK WITH POPULATIONS FOR VARIOUS INDEX TYPES:
C
C
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
C
C POSITIONS IN CORE WHERE THE GATHER-SCATTER VECTORS ARE LOCATED.
C
C        LEFT INDEX    - IRREP
C        CENTRAL INDEX - SYMMETRY TYPE
C        RIGHT INDEX   - ORBITAL TYPE (AI, IA, ETC)
C
      COMMON /FSIOLOC/ IFSPOS(8,22,8),ISCRLC, IFSPOSZ(8,22,8)
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
      DATA TYPE /'VRT','VRT' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','OCC' , 'VRT','VRT' , 'VRT','VRT' ,
     &           'OCC','OCC' , 'OCC','OCC' , 'VRT','OCC' ,
     &           'VRT','OCC' , 'VRT','OCC' , 'VRT','OCC' ,
     &           'VRT','VRT' , 'OCC','OCC' , 'VRT','VRT' ,
     &           'OCC','VRT' , 'OCC','VRT' , 'OCC','VRT' ,
     &           'VRT','VRT' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','OCC'/
      DATA SPIN / 1,1 , 2,2 , 1,1 , 2,2 , 1,1 , 2,2 , 1,1 , 2,2 , 
     &            1,1 , 2,2 ,
     &            1,2 , 2,1 , 1,2 , 1,2 , 1,2 , 1,1 , 2,2 , 1,2 ,
     &            1,1 , 2,2 ,
     &            1,1 , 2,2/
      DATA PACK /'UPLT','UPLT','UPLT','UPLT','UPLE',
     &           'UPLE','UPLE','UPLE','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL',
     &           'FULL','FULL'/
C SG 2/96
C      Data FSDsLbl /'AN', 'NA', 'AA', 'IN', 'NI', 'II', 'AI', 'IA'/
      FSDsLbl(1) = 'AN'
      FSDsLbl(2) = 'NA'
      FSDsLbl(3) = 'AA'
      FSDsLbl(4) = 'IN'
      FSDsLbl(5) = 'NI'
      FSDsLbl(6) = 'II'
      FSDsLbl(7) = 'AI'
      FSDsLbl(8) = 'IA'
C
C     SET TIOFull. If it is set to 1, the lists are 'right' after 
C             FSPUT and FSPUTT1. Otherwise only the active part
C	      can be used.
      TIOFull=TFull
C
C FILL VECTOR OF LENGTH NBAS WITH "1" 
C
      CALL IZERO(IFSPOS,8*22*8)
      CALL IZERO(IFSPOSZ,8*22*8)
      CALL IZERO(IVECTOR,MXBAS*6)
      DO 5 ISPIN=1,2
         DO 10 ITYPE=1,2
            DO 20 IPOS=1,NBAS
               IVECTOR(IPOS,ITYPE,ISPIN)=1
 20         CONTINUE
 10      CONTINUE
 5    CONTINUE
C
C     NOW ZERO OUT POSITIONS WHICH CORRESPOND TO ACTIVE ORBITALS (ITYPE=2)
C     AND POSITIONS WHICH CORRESPOND TO INACTIVE ORBITALS (ITYPE=3)
C
      DO 25 ISPIN=1,2
         DO 28 I=1,NACTIVE(ISPIN)
            IPOS=IACTIVE(I,ISPIN) 
            IVECTOR(IPOS,2,ISPIN)=0
            IVECTOR(IPOS,3,ISPIN)=1
 28      CONTINUE
 25   CONTINUE
C
      CALL FSINIPOP(IVECTOR,MXBAS,IUHF)
C
C     For ISYTYP 1-8, we need room for the "zero vector" accompanying 
C     the gather vector.  This is relatively small, and the space
C     required can be accurately and easily determined in advance:
C
C     AN + IN = NN, NA + NI = NN, AA + AI + IA + II = NN, or 3NN in total.
C
C     ISYTYP  Type    Expanded size
C        1    V <  V      V x V
C        2    v <  v      v x v
C        3    O <  O      O x O
C        4    o <  o      o x o
C        5    V <= V      V x V
C        6    v <= v      v x v
C        7    O <= O      O x O
C        8    o <= o      o x o
C     TOTAL:       2 ( VxV + vxv + OxO + oxo )
C
      IZOff = I0
C
      IZNeed = 3 * 2 * ( NVrtO(1)**2 + NVrtO(2)**2 + NOccO(1)**2 +
     $   NOccO(2)**2 )
      IOFF = I0 + IZNeed
cSSS      Write (6, *) '@FSIODRV-I, Allocated ', IZNeed, ' words for ',
cSSS     $   '''zero vectors''.'
C
C     CALCULATE THE LENGTH OF THE I/O BUFFER. ITS SPACE IS ALLOCATED
C     LATER IN THIS ROUTINE
C     THE ACTUAL SIZE IS GREATER BY ONE IN ORDER TO ALLOCATE SPACE
C     FOR THE 'DIAGONAL' ELEMENTS
C
      DSZMAX=0
      DO 100 I=1,22
       DO 101 J=1,8
        DSZMAX=MAX(DSZMAX,IRPDPD(J,I))
101    CONTINUE
100   CONTINUE
      DSZMAX=DSZMAX+1
C
C NOW GO THROUGH THE DIFFERENT DISTRIBUTION TYPES ONE BY ONE
C
      DO 30 ITYPE=1,22
       IF(IDO(ITYPE).EQ.1)THEN
        LFTTYP=TYPE(1,ITYPE)
        LFTSPN=SPIN(1,ITYPE) 
        RHTTYP=TYPE(2,ITYPE)
        RHTSPN=SPIN(2,ITYPE)
        TYPEX=PACK(ITYPE)
        IF(LFTTYP.EQ.'OCC')CALL ICOPY(8,POP(1,LFTSPN),1,POPLFT,1)
        IF(LFTTYP.EQ.'VRT')CALL ICOPY(8,VRT(1,LFTSPN),1,POPLFT,1)
        IF(RHTTYP.EQ.'OCC')CALL ICOPY(8,POP(1,RHTSPN),1,POPRHT,1)
        IF(RHTTYP.EQ.'VRT')CALL ICOPY(8,VRT(1,RHTSPN),1,POPRHT,1)
C
C LEFT INDEX ACTIVE, RIGHT INDEX FULL
C
cSSS        Write (6, 9999) 1, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,3,LFTSPN),IVECTOR(1,1,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDAN(1,ITYPE),
     &             IFSPOS(1,ITYPE,1),IOFF,
     $     IFSPOSZ(1, ITYPE, 1), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX FULL  , RIGHT INDEX ACTIVE
C
cSSS        Write (6, 9999) 2, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,1,LFTSPN),IVECTOR(1,3,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDNA(1,ITYPE),
     &             IFSPOS(1,ITYPE,2),IOFF,
     $     IFSPOSZ(1, ITYPE, 2), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX ACTIVE  , RIGHT INDEX ACTIVE
C
cSSS        Write (6, 9999) 3, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,3,LFTSPN),IVECTOR(1,3,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDAA(1,ITYPE),
     &             IFSPOS(1,ITYPE,3),IOFF,
     $     IFSPOSZ(1, ITYPE, 3), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX INACTIVE, RIGHT INDEX FULL
C
cSSS        Write (6, 9999) 4, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,2,LFTSPN),IVECTOR(1,1,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDIN(1,ITYPE),
     &             IFSPOS(1,ITYPE,4),IOFF,
     $     IFSPOSZ(1, ITYPE, 4), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX FULL , RIGHT INDEX INACTIVE
C
cSSS        Write (6, 9999) 5, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,1,LFTSPN),IVECTOR(1,2,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDNI(1,ITYPE),
     &             IFSPOS(1,ITYPE,5),IOFF,
     $     IFSPOSZ(1, ITYPE, 5), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX INACTIVE, RIGHT INDEX INACTIVE
C
cSSS        Write (6, 9999) 6, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,2,LFTSPN),IVECTOR(1,2,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDII(1,ITYPE),
     &             IFSPOS(1,ITYPE,6),IOFF,
     $     IFSPOSZ(1, ITYPE, 6), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX ACTIVE, RIGHT INDEX INACTIVE
C
cSSS        Write (6, 9999) 7, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,3,LFTSPN),IVECTOR(1,2,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDAI(1,ITYPE),
     &             IFSPOS(1,ITYPE,7),IOFF,
     $     IFSPOSZ(1, ITYPE, 7), ICore(IZOff), Izoff,Dszmax)
C
C LEFT INDEX INACTIVE, RIGHT INDEX ACTIVE
C
cSSS        Write (6, 9999) 8, IType, LftTyp, LftSpn, RhtTyp, RhtSpn, TypeX
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,2,LFTSPN),IVECTOR(1,3,RHTSPN),NBAS,
     &             ICORE(IOFF),FSDPDIA(1,ITYPE),
     &             IFSPOS(1,ITYPE,8),IOFF,
     $     IFSPOSZ(1, ITYPE, 8), ICore(IZOff), Izoff,Dszmax)
       ENDIF
30    CONTINUE
 9999 Format(1X, '@FSIODRV-I, distribution type ', I1, ' list type ',
     $   I2, ': ', A, '(', I1, '), ', A, '(', I1, ') ', A)
C
      WRITE(6,1000) '''zero''', IZOFF-I0
      WRITE(6,1000) 'I/O', IOff-IZOFF
1000  FORMAT(T3,'@FSIODRV-I, ', A, ' vectors require ',I8,' words.')
C
C NOW ALLOCATE SPACE FOR THE LARGEST DISTRIBUTION.  THIS IS USED
C AS A SCRATCH AREA IN FSGET. 
C
      IOFF=IOFF+1-iand(IOFF,1)
      ISCRLC=IOFF
      MAXCOR=MAXCOR-(IOFF-I0)-DSZMAX*IINTFP
      I0=IOFF+DSZMAX*IINTFP
      TIOBufSz = DSzMax
C
      Write (6, 1020) DSzMax * IIntFP
 1020 Format(T3, '@FSIODRV-I, I/O buffer requires ', I8, ' words.')
C
C     Find out how much to allocate for "temporary" I/O lists.
C     Putting this in a subroutine reduces clutter here and
C     localizes things that have to change for "temorary" I/O.
C
      Call AllocTIO(MaxCor - I0, TIOSiz)
      Write(6, 1010) TIOSiz
 1010 Format(T3, '@FSIODRV-I, Temporary I/O vectors allocated ', I8,
     $   ' words.')
      TIOBas = I0
      I0 = I0 + TIOSiz
      MAXCOR=MAXCOR-TIOSiz
C
      RETURN
      END
