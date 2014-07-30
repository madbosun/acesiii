      SUBROUTINE FSMAP(XTYP,YTYP,XSPIN,YSPIN,TYPE,POPX,POPY,IDROPX,
     &                 IDROPY,NBAS,IMAP,DPDVEC,IFSPOS,
     &                 IOFF, IZPos, ZList, IZOff,DSZMAX)
C
C THIS ROUTINE CONSTRUCTS A GATHER-SCATTER VECTOR WHICH RELATES
C POSITIONS IN A SYMMETRY PACKED PQ DISTRIBUTION TO THOSE IN THE
C SAME DISTRIBUTION WHERE CERTAIN ELEMENTS OF P OR Q ARE EXCLUDED.
C
C THIS VECTORS PRODUCED BY THIS ROUTINE ARE USEFUL FOR READING QUANTITIES
C INVOLVING ONLY ACTIVE OR INACTIVE INDICES FROM THE STANDARD LISTS.
C
C The last three arguments are used only for TYPE = 'UPLT' or 'UPLE'
C (unpack from X < Y  or X <= Y storage) 
C
C     IZPos   Analog to IFSPos for ZList. 
C     ZList   Integer multiplier for each element of IMap.
C     IZOff   Analog to IOff for tracking the position of the ZList
C     DSZMAX  LENGTH OF THE FSIO BUFFER
C
C NOTE: DPDVEC applies equally well to IMap and ZList.
C
CEND
      IMPLICIT INTEGER (A-Z)
C
      CHARACTER*4 TYPE
      CHARACTER*3 XTYP,YTYP
      DIMENSION POPX(8),POPY(8),IDROPX(NBAS),IDROPY(NBAS),IMAP(1)
      DIMENSION IOFFX(8),IOFFY(8),DPDVEC(8),IFSPOS(8)
      Integer IZPos(8), ZList(*), IZOff
C
C     Local Variables
C
      Integer Ptrs(8,2), Count, I
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /INFO  / NOCCO(2),NVRTO(2)
C
      IPART=0
C
C CALCULATE OFFSETS FOR STARTS OF IRREPS
C
      IF(XTYP.EQ.'VRT')THEN
       IOFFX(1)=NOCCO(XSPIN)
      ELSE
       IOFFX(1)=0
      ENDIF
      IF(YTYP.EQ.'VRT')THEN
       IOFFY(1)=NOCCO(YSPIN)
      ELSE
       IOFFY(1)=0
      ENDIF
      DO 5 IRREP=2,NIRREP
       IOFFX(IRREP)=IOFFX(IRREP-1)+POPX(IRREP-1)
       IOFFY(IRREP)=IOFFY(IRREP-1)+POPY(IRREP-1) 
5     CONTINUE
      IABSX=0
      IABSY=0
C
      IF(TYPE.EQ.'FULL')THEN
C
C X,Y (FULL) DISTRIBUTION.  
C
       DO 10 IRREPXY=1,NIRREP 
        INUMIRR=0
        IFULL0=0
        DO 20 IRREPY=1,NIRREP
         IRREPX=DIRPRD(IRREPY,IRREPXY)
         NUMX=POPX(IRREPX)
         NUMY=POPY(IRREPY)
         DO 30 INDEXY=1,NUMY
          DO 40 INDEXX=1,NUMX
           IABSX=IOFFX(IRREPX)+INDEXX
           IABSY=IOFFY(IRREPY)+INDEXY
           IFULL0=IFULL0+1
           IF(IDROPX(IABSX)*IDROPY(IABSY).NE.0)THEN
            IPART=IPART+1
            INUMIRR=INUMIRR+1
            IMAP(IPART)=IFULL0
           ENDIF
40        CONTINUE
30       CONTINUE 
20      CONTINUE
        DPDVEC(IRREPXY)=INUMIRR
        IFSPOS(IRREPXY)=IOFF
        IOFF=IOFF+INUMIRR
10     CONTINUE
C
      ELSEIF(TYPE.EQ.'PACK'.OR.TYPE.EQ.'PCK2')THEN
C
C X<Y OR X<=Y (PACKED) DISTRIBUTION 
C
       IF(TYPE.EQ.'PACK')THEN
        IOFFR=1
       ELSE
        IOFFR=0
       ENDIF
       DO 110 IRREPXY=1,NIRREP
        INUMIRR=0
        IFULL0 =0
        DO 120 IRREPY=1,NIRREP
         IRREPX=DIRPRD(IRREPY,IRREPXY)
C
         IF(IRREPX.EQ.IRREPY)THEN
C
          NUMX=POPX(IRREPX)
          NUMY=POPY(IRREPY)
          DO 130 INDEXY=1,NUMX
           DO 140 INDEXX=1,INDEXY-IOFFR
            IABSX=IOFFX(IRREPX)+INDEXX
            IABSY=IOFFY(IRREPY)+INDEXY
            IFULL0=IFULL0+1
            IF(IDROPX(IABSX)*IDROPY(IABSY).NE.0)THEN
             IPART=IPART+1
             INUMIRR=INUMIRR+1
             IMAP(IPART)=IFULL0
            ENDIF
140        CONTINUE
130       CONTINUE
C
         ELSEIF(IRREPX.LT.IRREPY)THEN
C          
          NUMX=POPX(IRREPX)
          NUMY=POPY(IRREPY)
          DO 135 INDEXY=1,NUMY
           DO 145 INDEXX=1,NUMX
            IABSX=IOFFX(IRREPX)+INDEXX
            IABSY=IOFFY(IRREPY)+INDEXY
            IFULL0=IFULL0+1
            IF(IDROPX(IABSX)*IDROPY(IABSY).NE.0)THEN
             IPART=IPART+1
             INUMIRR=INUMIRR+1
             IMAP(IPART)=IFULL0
            ENDIF
145        CONTINUE
135       CONTINUE 
C
         ENDIF
C
120     CONTINUE
        DPDVEC(IRREPXY)=INUMIRR
        IFSPOS(IRREPXY)=IOFF
        IOFF=IOFF+INUMIRR
110    CONTINUE
C
      ELSEIF(TYPE.EQ.'UPLE'.OR.TYPE.EQ.'UPLT')THEN
C
C        X<=Y OR X<Y (PACKED) DISTRIBUTION which is to be expanded
C        into a rectangle.  The difference between this and 'PACK' &
C        'PCK2' above is that we have to cover the entire array,
C        to insure that things expand properly, but the origin
C        position in IMap must be within the X<Y or X<=Y framework.
C
C        Things to be very careful of in comparison with regular
C        PACK or PCK2 situations:
C
C        o For X<=Y lists, diagonals (X=Y) are not stored.  They are
C          of course required in the unpacked list.
C        o Integrals, amplitudes, etc. are antisymmetric, so unpacking
C          requires a sign change on some of elements compared to
C          what they were in the packed list.
C
C        These points are solved by creating a companion array to IMap
C        which contains an integer multipliers for each element in
C        in the unpacked list.  Elements may be 1, 0, -1 for UPLT or
C        1, -1 for UPLE.  All offsets for the IMap array apply equally
C        to the Zlist array.
C
C        o The IrrepX .gt. IrrepY blocks of the unpacked list are
C          the transpose of the corresponding IrrepX <--> IrrepY
C          block of the packed list.  They must be given IMap entries
C          before the IFull0 value of the beginning of the actual
C          block in the packed list is known because each element of
C          IMap corresponds directly to an element in the unpacked
C          result, so proper ordering must be maintained.
C
C        This is solved by filling the IrrepX .gt. IrrepY blocks with
C        the position relative to the beginning of the actual block
C        and saving the position and length of this block in the
C        IrrepX element of the Ptrs array.
C
C        When we process a block which is in the packed list (IrrepX
C        .lt. IrrepY), we lookup the IrrepY value in the Ptrs array
C        to find out where the transpose block is which requires
C        the value of IFull0 to be added.
C
         IF(TYPE.EQ.'UPLT')THEN
            IOFFR=1
         ELSE
            IOFFR=0
         ENDIF
         DO 210 IRREPXY=1,NIRREP
            INUMIRR=0
            NZero = 0
            IFULL0 =0
C
C           The pointer array is reused for each IRREPXY
C
            Call IZero(Ptrs, 8*2)
C
            DO 220 IRREPY=1,NIRREP
               IRREPX=DIRPRD(IRREPY,IRREPXY)
C
               IF(IRREPX.EQ.IRREPY)THEN
C
                  NUMX=POPX(IRREPX)
                  NUMY=POPY(IRREPY)
                  DO 230 INDEXY=1,NUMY
                     DO 240 INDEXX=1,INDEXY - IOffR
                        IABSX=IOFFX(IRREPX)+INDEXX
                        IABSY=IOFFY(IRREPY)+INDEXY
                        IFULL0=IFULL0+1
                        IF(IDROPX(IABSX)*IDROPY(IABSY).NE.0)THEN
                           IPART=IPART+1
                           INUMIRR=INUMIRR+1
                           IMAP(IPART)=IFULL0
                           Zlist(IPart) = 1
                        ENDIF
 240                 CONTINUE
C
                     IABSX=IOFFX(IRREPX)+INDEXX+1
                     IABSY=IOFFY(IRREPY)+INDEXY
                     If (Type .eq. 'UPLT' .AND. 
     $                  (IDropX(IAbsY)*IDropY(IAbsY) .ne. 0) ) then
                        IPart = IPart + 1
                        INumIrr = INumIrr + 1
c&C
c&C                       In case this is the (1,1) element, make sure
c&C                       it is at least 1
c&C
c&                        IMap(IPart) = Max(1, IFull0)
                        IMap(IPart) = DSZMAX
C
C                       This one will have to be explicitly zeroed in
c                       the final gathered array.  DO NOT increment
c                       IFull0 here!
C
                        ZList(IPart) = 0
                     EndIf
C
C                    These are the transpose of the upper triangle
C                    so a sign change is required.
C
C                    Increment to IFull1 depends on whether list is
C                    UPLE or UPLT.
C
                     IFull1 = IFull0 + IndexY
                     DO 250 INDEXX=INDEXY+1, NUMX
                        IABSX=IOFFX(IRREPX)+INDEXX
                        IABSY=IOFFY(IRREPY)+INDEXY
cSSS                        IFULL1=IFULL0 + IndexX - 1
                        IF(IDROPX(IABSX)*IDROPY(IABSY).NE.0)THEN
                           IPART=IPART+1
                           INUMIRR=INUMIRR+1
                           IMAP(IPART)=IFULL1
                           ZList(IPart) = -1
                        ENDIF
                        IFull1 = IFull1 + IndexX - IOffR
 250                 CONTINUE
 230              CONTINUE
C                 
               ElseIf (IrrepX .gt. IrrepY) then
C              
C                 This block is the transpose of the IrrepY,IrrepX block
C                 which will satisfy the conditional below.  At the
C                 moment, we don't know where the original block will
C                 begin, so just put offsets into IMap now and setup
C                 Ptrs so that when we hit the original block, we can
C                 add the appropriate offset to these values.
C
                  Count = 0
                  Ptrs(IrrepX, 1) = IPart + 1
C
                  NumX = PopX(IrrepX)
                  NumY = PopY(IrrepY)
                  Do 232 IndexY = 1, NumY
                     Do 242 IndexX = 1, NumX
                        IAbsX = IOffX(IrrepX) + IndexX
                        IAbsY = IOffY(IrrepY) + IndexY
C
C                       Remember: original position is in the
C                       _transpose_ of this block w.r.t. the
C                       current NumX, NumY.  The formula may
C                       look funny, but its correct -- try it!
C
                        If (IDropX(IAbsX) * IDropY(IAbsY) .ne. 0) then
                           Count = Count + 1
                           IMap(IPart + Count) =
     $                        (IndexX-1) * NumY + IndexY
                           ZList(IPart + Count) = -1
                        EndIf
C
 242                  Continue
 232               Continue
C
                   Ptrs(IrrepX, 2) = Count
                   IPart = IPart + Count
                   INumIrr = INumIrr + Count
C                 
               ELSEIF(IRREPX.LT.IRREPY)THEN
C          
C                 Check for a block in the lower half of the matrix
C                 which needs to know the origin of this block in
C                 the original triangle.
C
                  If (Ptrs(IrrepY, 2) .ne. 0) then
                     Do 234 I = 0, Ptrs(IrrepY, 2) - 1
                        J = Ptrs(IrrepY, 1) + I
                        IMap(J) = IMap(J) + IFull0
 234                 Continue
                  EndIf
C
C                 Now continue with this block
C
                  NUMX=POPX(IRREPX)
                  NUMY=POPY(IRREPY)
                  DO 235 INDEXY=1,NUMY
                     DO 245 INDEXX=1,NUMX
                        IABSX=IOFFX(IRREPX)+INDEXX
                        IABSY=IOFFY(IRREPY)+INDEXY
                        IFULL0=IFULL0+1
                        IF(IDROPX(IABSX)*IDROPY(IABSY).NE.0)THEN
                           IPART=IPART+1
                           INUMIRR=INUMIRR+1
                           IMAP(IPART)=IFULL0
                           ZList(IPart) = 1
                        ENDIF
 245                 CONTINUE
 235              CONTINUE 
C                 
               ENDIF
C
 220        CONTINUE
C
            DPDVEC(IRREPXY)=INUMIRR
            IFSPOS(IRREPXY)=IOFF
            IOFF=IOFF+INUMIRR
C           
            If (Type .eq. 'UPLT' .OR. Type .eq. 'UPLE') then
               IZPos(IrrepXY) = IZOff
               IZOff = IZOff + INumIrr
            EndIf
 210     CONTINUE
C
      ENDIF
C
cSSS      If (Type .eq. 'UPLT' .OR. Type .eq. 'UPLE') then
cSSS         Write (6, 9005) (i, IFSPOS(i), IZPos(i), DPDVEC(i),
cSSS     $      i = 1, NIrrep)
cSSS      Else
cSSS         Write (6, 9006) (i, IFSPOS(i), DPDVEC(i), i = 1, NIrrep)
cSSS      EndIf
cSSS      If (Type .eq. 'UPLT' .OR. Type .eq. 'UPLE') then
cSSS         Write (6, 9011) (i, Zlist(i), IMap(I), ZList(i+1), IMap(i+1),
cSSS     $      ZList(I+2), IMap(I+2), ZList(I+3), IMap(I+3),
cSSS     $      ZList(I+4), IMap(i+4), I = 1, IPart, 5)
cSSS      Else
cSSS         Write (6, 9010) (i, IMap(I), IMap(i+1), IMap(I+2), IMap(I+3),
cSSS     $   IMap(i+4), I = 1, IPart, 5)
cSSS      EndIf
cSSS 9005 Format(10X, 'IFSPOS', 2X, 'IZPOS ', 2X, 'DPDVEC', /
cSSS     $   (5X, I3, 2X, I6, 2I8))
cSSS 9006 Format(10X, 'IFSPOS', 2X, 'DPDVEC', /
cSSS     $   (5X, I3, 2X, I6, I8))
cSSS 9010 Format(/ (5X, I3, ':', 2X, 5I6) )
cSSS 9011 Format(/ (5X, I3, ':', 2X, 5(I2, '*', I4, 4X) ) ) 
C
      RETURN
      END
