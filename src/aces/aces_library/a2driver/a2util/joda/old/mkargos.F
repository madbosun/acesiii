      SUBROUTINE MKARGOS(AGEN,NGENG,PGRP,NTATOM,ZSYM,BASNAM,NORD,COORD,
     +         TITLE,IATNUM,ITAPE,NREALAT,IMAP)
C
C     WRITES ARGOS INPUT
C
C     WRITEN BY P.G.SZALAY (JUL 92)
C     BASED ON THE AUTBAS PROGRAM BY STANISLAV BISKUPIC AND OTHERS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

csb 1/97 Hold long basis set names
#include "mxatms.par"
      integer baslen
      parameter (baslen=80)
      character*(baslen) ibas(mxatms)
      
      CHARACTER*(*) PGRP,ZSYM(*),BASNAM(*),TITLE
      CHARACTER*60 IATM
      CHARACTER*1 ACHAR,AGEN
      DIMENSION AGEN(3,3)
      DIMENSION COORD(*),NORD(*),IATNUM(*),IMAP(NREALAT)
      DIMENSION NSYT(mxatms),IATM(mxatms)
      CHARACTER*4 ICH1,ICH2,ICH3,ICH4,ICH5,ICH6,ICH7,ICH8,IBF
      DIMENSION IBF(8)
      DIMENSION ICH2(2),ICH3(2),ICH4(2),ICH5(4),ICH6(4),ICH7(4),ICH8(8)
      DATA ICH2/'A''','A"'/,ICH3/'A','B'/,ICH4/'AG','BU'/
      DATA ICH5/'A','B3','B2','B1'/,ICH6/'A1','B1','B2','A2'/
      DATA ICH7/'AG','AU','BG','BU'/,ICH1/'A'/
      DATA ICH8/'AG','B3U','B2U','B1G','B1U','B2G','B3G','AU'/
      INTEGER SHELLANG,SHELLLOC,SHELLSIZ,SHELLPRM
      COMMON /MOLSTR1/ SHELLANG(100),SHELLLOC(100),SHELLSIZ(100),
     &                 SHELLPRM(100),BASISEXP(10000),BASISCNT(10000)
      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IATOM
      COMMON /FLAGS/IFLAGS(100),IFLAGS2(500)
      COMMON/WORK/  IDUM(10000),INFO(100,3)
      DIMENSION ICHHX(8,3),NBFS(8)
      DOUBLE PRECISION RR(3,200),ATNM(150)
      DIMENSION IATN(150),JAOR(5)
      COMMON/JJJ/ JJQQ(8),JJPP(8)
      DIMENSION IBING(3,3),ICHAT(8,8),ICHGEN(8,3),IMAGE(200,3)
      COMMON /IB/IBING,ICHAT,ICHGEN,IMAGE,NEQA
      DIMENSION ICHAS(3,1),ICHAP(3,3),ICHAD(3,6),ICHAF(3,10),ICHAG(3,15)
      DIMENSION ACHAR(3),ISOM(6400),ICHSO(100),NEQ(200),LOGO(3)
      DIMENSION NSTV(3,7),NNSTV(8),ICCN(30)
      DIMENSION RNUC(3,8),RSCR(3)
      LOGICAL ISORB,IPORB,IDORB,IFORB,IGORB,LOGO
      COMMON/JODALOG/ ISORB,IPORB,IDORB,IFORB,IGORB
      DIMENSION JCH(8)
C
      COMMON/IRD/ IRD(5),IFL1(100)
C
      DATA JCH/1,2,2,2,4,4,4,8/
      DIMENSION IPRMS(1,1)
      DATA IPRMS/1/
      DIMENSION IPRMP(3,3)
      DATA IPRMP/1,0,0,
     &           0,1,0,
     &           0,0,1/
      DIMENSION IPRMD(6,5)
      DATA IPRMD/-1,-1, 4, 0, 0, 0,
     &            0, 0, 0, 0, 1, 0,
     &            0, 0, 0, 0, 0, 1,
     &            0, 0, 0, 1, 0, 0,
     &            1,-1, 0, 0, 0, 0/
      DIMENSION IPRMDC(6,6)
      DATA IPRMDC/-1,-1, 4, 0, 0, 0,
     &             0, 0, 0, 0, 1, 0,
     &             0, 0, 0, 0, 0, 1,
     &             0, 0, 0, 1, 0, 0,
     &             1,-1, 0, 0, 0, 0,
     &             1, 1, 1, 0, 0, 0/
      DIMENSION IPRMF(10,7)
      DATA IPRMF/ 0, 0, 4, 0,-9, 0,-9, 0, 0, 0,
     &           -1, 0, 0, 0, 0,-1, 0,16, 0, 0,
     &            0,-1, 0,-1, 0, 0, 0, 0,16, 0,
     &            0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &            0, 0, 0, 0, 1, 0,-1, 0, 0, 0,
     &            1, 0, 0, 0, 0,-9, 0, 0, 0, 0,
     &            0, 1, 0,-9, 0, 0, 0, 0, 0, 0/
      DIMENSION IPRMFC(10,10)
      DATA IPRMFC/ 0, 0, 4, 0,-9, 0,-9, 0, 0, 0,
     &            -1, 0, 0, 0, 0,-1, 0,16, 0, 0,
     &             0,-1, 0,-1, 0, 0, 0, 0,16, 0,
     &             0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &             0, 0, 0, 0, 1, 0,-1, 0, 0, 0,
     &             1, 0, 0, 0, 0,-9, 0, 0, 0, 0,
     &             0, 1, 0,-9, 0, 0, 0, 0, 0, 0,
     &             1, 0, 0, 0, 0, 1, 0, 1, 0, 0,
     &             0, 1, 0, 1, 0, 0, 0, 0, 1, 0,
     &             0, 0, 1, 0, 1, 0, 1, 0, 0, 0/
      DIMENSION IPRMG(15,9)
      DATA IPRMG/9, 9,64, 0, 0, 0, 0, 0, 0, 36,-576,-576, 0, 0, 0,
     &           0, 0, 0, 0,-9, 0, 0,16, 0,  0,   0,   0, 0,-9, 0,
     &           0, 0, 0, 0, 0, 0,-9, 0,16,  0,   0,   0,-9, 0, 0,
     &           0, 0, 0,-1, 0,-1, 0, 0, 0,  0,   0,   0, 0, 0,36,
     &           1,-1, 0, 0, 0, 0, 0, 0, 0,  0, -36,  36, 0, 0, 0,
     &           0, 0, 0, 0, 1, 0, 0, 0, 0,  0,   0,   0, 0,-9, 0,
     &           0, 0, 0, 0, 0, 0, 1, 0, 0,  0,   0,   0,-9, 0, 0,
     &           0, 0, 0, 1, 0,-1, 0, 0, 0,  0,   0,   0, 0, 0, 0,
     &           1, 1, 0, 0, 0, 0, 0, 0, 0,-36,   0,   0, 0, 0, 0/
      DIMENSION IPRMGC(15,15)
      DATA IPRMGC/ 1, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,
     &             0, 1, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,
     &             0, 0, 1, 0, 0, 0,0,0,0,0,0,0,0,0,0,
     &             0, 0, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0,
     &             0, 0, 0, 0, 1, 0,0,0,0,0,0,0,0,0,0,
     &             0, 0, 0, 0, 0, 1,0,0,0,0,0,0,0,0,0,
     &             0, 0, 0, 0, 0, 0,1,0,0,0,0,0,0,0,0,
     &             0, 0, 0, 0, 0, 0,0,1,0,0,0,0,0,0,0,
     &             0, 0, 0, 0, 0, 0,0,0,1,0,0,0,0,0,0,
     &             0, 0, 0, 0, 0, 0,0,0,0,1,0,0,0,0,0,
     &             0, 0, 0, 0, 0, 0,0,0,0,0,1,0,0,0,0,
     &             0, 0, 0, 0, 0, 0,0,0,0,0,0,1,0,0,0,
     &             0, 0, 0, 0, 0, 0,0,0,0,0,0,0,1,0,0,
     &             0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,1,0,
     &             0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,1/
C
      DATA ICHAS,ICHAP,ICHAD/3*1,-1,3*1,-1,3*1,-1,9*1,-1,-1,1,-1,
     11,-1,1,-1,-1/
      DATA ICHAF/-1,3*1,-1,3*1,-1,1,-1,3*1,-1,-1,4*1,-1,-1,3*1,-1,1,
     13*-1/
      DATA ICHAG/9*1,2*-1,1,-1,1,3*-1,2*1,3*-1,1,-1,1,2*-1,10*1,
     &           3*-1,1,3*-1,1/
      DATA NSTV/4,2,3,6,2,5,7,3,5,7,6,4,8,4,5,8,6,3,8,7,2/
      DATA ACHAR/'X','Y','Z'/,NNSTV/4*0,3*1,7/
      DIMENSION ITABH(4,4),ICHGH(4,3)
      DATA ITABH/4*1,1,-1,1,-1,1,1,-1,-1,1,-1,-1,1/
      DATA ICHGH/1,-1,-1,1,1,-1,1,-1,1,-1,-1,1/
      DATA ICHHX/1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,1,1,-1,-1,
     1  1,1,1,1,-1,-1,-1,-1/
C
      CALL IZERO(IMAP,NREALAT)
      DO 1948 K=1,150
 1948 IATN(K)=0
      DO 1949 I=1,8
      DO 1949 J=1,3
 1949 ICHGEN(I,J)=ICHHX(I,J)
      IAORDS=0
      ICCX=0
      DO 205 I=1,20
      DO 205 J=1,3
  205 IMAGE(I,J)=0
C
      IF(NGENG.EQ.0) GO TO 175
C
C     CONVERT SYMBOLIC GENERATORS AGEN TO BINARY FORM IBING
C     AGEN(I) = X**IBING(1,I) Y**IBING(2,I) Z**IBING(3,I)
C
      DO 151 IJ=1,3
      DO 151 JI=1,3
  151 IBING(IJ,JI)=0
      DO 170 I=1,NGENG
      DO 150 L=1,3
  150 IBING(L,I)=0
      DO 170 J=1,3
      DO 160 K=1,3
      IF(AGEN(J,I).EQ.ACHAR(K)) IBING(K,I)=1
  160 CONTINUE
  170 CONTINUE
  175 CONTINUE
      IGR=INTGRP(PGRP)
      JJCH=JCH(IGR)
      GO TO(3451,3452,3453,3454,3455,3456,3457,3458),IGR
 3451 CONTINUE
      IBF(1)=ICH1
      GO TO 3459
 3452 CONTINUE
      DO 182 KL=1,JJCH
  182 IBF(KL)=ICH2(KL)
      GO TO 3459
 3453 CONTINUE
      DO 183 KL=1,JJCH
  183 IBF(KL)=ICH3(KL)
      GO TO 3459
 3454 CONTINUE
      DO 184 KL=1,JJCH
  184 IBF(KL)=ICH4(KL)
      GO TO 3459
 3455 CONTINUE
      DO 185 KL=1,JJCH
  185 IBF(KL)=ICH5(KL)
      GO TO 3459
 3456 CONTINUE
      DO 186 KL=1,JJCH
  186 IBF(KL)=ICH6(KL)
      GO TO 3459
 3457 CONTINUE
      DO 187 KL=1,JJCH
  187 IBF(KL)=ICH7(KL)
      GO TO 3459
 3458 CONTINUE
      DO 188 KL=1,JJCH
  188 IBF(KL)=ICH8(KL)
 3459 CONTINUE
      IF(IGR.NE.2.AND.IGR.NE.7) GO TO 101
      DO 102 I=1,4
      DO 102 J=1,3
  102 ICHGEN(I,J)=ICHGH(I,J)
  101 CONTINUE
      DO 3219 I9=1,8
      JJPP(I9)=I9
 3219 JJQQ(I9)=I9
C REORDER OF IRREPS. IN PRINCIPLE WE COULD DO IT USING OCCUPATIONS
C&      IF(IGR.EQ.1) GO TO 181
C&      CALL REORD(IBF,JJCH,ICHGEN,NSTV,NNSTV(IGR))
C&  181 CONTINUE
      IF(IGR.NE.2.AND.IGR.NE.7) GO TO 103
      DO 104 I=1,4
      DO 104 J=1,4
  104 ICHAT(I,J)=ITABH(I,J)
  103 DO 130 N=1,8
      II=0
      DO 120 I=1,2
      DO 110 J=1,2
      DO 100 K=1,2
      II=II+1
  100 ICHAT(N,II)=ICHGEN(N,1)**(K-1)*ICHGEN(N,2)**(J-1)*ICHGEN(N,3)
     1 **(I-1)
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
C
C SET DEFAULTS AND OFFSETS
C
      NGEN=0
      NTOTSHEL=0
      INS=0
      KR=0
      IK=0
      NATOMS=0
      ISST=0
      IOFFPRIM=1
      IOFFSHEL=1
C
      DO 300 ITATOM=1,NTATOM
       INS=INS+1
       NStart=1
       Do 319 J=1,Itatom-1
        NStart=NStart+Nord(MxAtms+J)
 319   Continue
       IBAS(INS)=BASNAM(NORD(NSTART))
       IATTYP=1
       IIATOM=NORD(NSTART)
       IATOM=IATNUM(IIATOM)
       IF(IATOM.GT.2)IATTYP=2
       IF(IATOM.GT.10)IATTYP=3 
       IATM(INS)=ZSYM(IIATOM)
       ATNM(INS)=FLOAT(IATOM)
C     
       ISTR=0
       DO 290 ISYT=1,INS-1
          IF(IBAS(ISYT).EQ.IBAS(INS)) THEN
             ICCN(INS)=ICCN(ISYT)
             ISTR=1
             GOTO 2901
          ENDIF
 290   CONTINUE
       ICCN(INS)=NTOTSHEL
 2901  CONTINUE
C   
       CALL RDBASPS(IBAS(INS),IFLAGS(64+IATTYP),ISTR,NSHELL)
       NSYT(INS)=NSHELL
C
       IBOTOM=3*NORD(NSTART)-2
       DO 291 ICOR=1,3
          RNUC(ICOR,1)=COORD(IBOTOM+ICOR-1)
 291   CONTINUE
          CALL ARGZMAP(RNUC(1,1),IMAPC,COORD,NREALAT)
          IF(IMAPC.NE.NORD(NSTART)) THEN
             WRITE(6,*)'@MKARG-E Internal error calling ARGZMAP'
             CALL ERREX
          ENDIF
       NEQA=1
       IF (NGENG.EQ.0) GO TO 202
       DO 200 I=1,NGENG
          DO 190 J=1,NEQA
             RD=0
             DO 180 K=1,3
                RNUC(K,NEQA+J)=RNUC(K,J)*(-1)**IBING(K,I)
                RD=RD+DABS(RNUC(K,1)-RNUC(K,NEQA+J))
 180         CONTINUE
             IF (RD.LT.1.D-10) GO TO 195
 190      CONTINUE
          NEQA=NEQA*2
 195      CONTINUE
 200   CONTINUE
 202   CONTINUE
C
       ISST=ISST+1
       NEQ(ISST)=NEQA
       NATOMS=NATOMS+NEQA
       DO 203 I=1,NEQA
          KR=KR+1
          DO 204 J=1,3
             RR(J,KR)=RNUC(J,I)
 204      CONTINUE
          CALL ARGZMAP(RR(1,KR),IMAP(KR),COORD,NREALAT)
 203   CONTINUE
C
C     DETERMINE IMAGING PATTERN OF NUCLEI UNDER GENERATORS
C
       IF (NGENG.EQ.0) GOTO 3906
       DO 240 I=1,NEQA
          IK=IK+1
          DO 250 IGEN =1,NGENG
             DO 210 K=1,3
                RSCR(K)=RNUC(K,I)*(-1)**IBING(K,IGEN)
 210         CONTINUE
             DO 230 J=1,NEQA
                RD=0
                DO 220 K=1,3
                   RD=RD+DABS(RNUC(K,J)-RSCR(K))
 220            CONTINUE
                IF (RD.LT.1.D-10) IMAGE(IK,IGEN)=J
 230         CONTINUE
 250      CONTINUE
 240   CONTINUE
 3906  CONTINUE
C
       IF(.NOT.ISORB) GO TO 401
       CALL MSOMAT(IPRMS,1,1,ISOM,NEQA,NEQA,ICHAS,ICHSO,NATOMS,NGENG)
       IAORDS=IAORDS+1
       CALL CHECK(IAORDS,ISOM,ICHSO,NEQA,NEQA,JAOR(1))
       IRDS=IRD(1)
       IF(IRDS.EQ.0) GO TO 8001
       DO 1661 I5=1,IRDS
          ICCX=ICCX+1
          IFL1(ICCX)=JAOR(1)
 1661  CONTINUE
 8001  CONTINUE
C
 401   IF(.NOT.IPORB) GO TO 402
       CALL MSOMAT(IPRMP,3,3,ISOM,3*NEQA,3*NEQA,ICHAP,ICHSO,NATOMS,
     &    NGENG)
       IAORDS=IAORDS+1
       CALL CHECK(IAORDS,ISOM,ICHSO,3*NEQA,3*NEQA,JAOR(2))
       IRDP=IRD(2)
       IF(IRDP.EQ.0) GO TO 8002
       DO 1662 I5=1,IRDP
          ICCX=ICCX+1
            IFL1(ICCX)=JAOR(2)
 1662  CONTINUE
 8002  CONTINUE
C
 402   IF(.NOT.IDORB) GO TO 403
       IF(IFLAGS(62).EQ.1) THEN
         CALL MSOMAT(IPRMD,6,5,ISOM,5*NEQA,6*NEQA,ICHAD,ICHSO,NATOMS,
     &       NGENG)
         IAORDS=IAORDS+1
         CALL CHECK(IAORDS,ISOM,ICHSO,5*NEQA,6*NEQA,JAOR(3))
      ELSE
         CALL MSOMAT(IPRMDC,6,6,ISOM,6*NEQA,6*NEQA,ICHAD,ICHSO,NATOMS,
     &       NGENG)
         IAORDS=IAORDS+1
         CALL CHECK(IAORDS,ISOM,ICHSO,6*NEQA,6*NEQA,JAOR(3))
      ENDIF
       IRDD=IRD(3)
       IF(IRDD.EQ.0) GO TO 8003
       DO 1663 I5=1,IRDD
          ICCX=ICCX+1
          IFL1(ICCX)=JAOR(3)
 1663  CONTINUE
 8003  CONTINUE
C
 403   IF(.NOT.IFORB) GO TO 404
       IF(IFLAGS(62).EQ.1) THEN
          CALL MSOMAT(IPRMF,10,7,ISOM,7*NEQA,10*NEQA,ICHAF,ICHSO,
     +       NATOMS,NGENG)
          IAORDS=IAORDS+1
          CALL CHECK(IAORDS,ISOM,ICHSO,7*NEQA,10*NEQA,JAOR(4))
       ELSE
          CALL MSOMAT(IPRMFC,10,10,ISOM,10*NEQA,10*NEQA,ICHAF,ICHSO,
     +       NATOMS,NGENG)
          IAORDS=IAORDS+1
          CALL CHECK(IAORDS,ISOM,ICHSO,10*NEQA,10*NEQA,JAOR(4))
       ENDIF
       IRDF=IRD(4)
       IF(IRDF.EQ.0) GO TO 8004
       DO 1664 I5=1,IRDF
          ICCX=ICCX+1
          IFL1(ICCX)=JAOR(4)
 1664  CONTINUE
 8004  CONTINUE
C
 404   IF(.NOT.IGORB) GO TO 405
       IF(IFLAGS(62).EQ.1) THEN
         CALL MSOMAT(IPRMG,15,9,ISOM,9*NEQA,15*NEQA,ICHAG,ICHSO,
     +       NATOMS,NGENG)
         IAORDS=IAORDS+1
         CALL CHECK(IAORDS,ISOM,ICHSO,9*NEQA,15*NEQA,JAOR(5))
       ELSE
         CALL MSOMAT(IPRMGC,15,15,ISOM,15*NEQA,15*NEQA,ICHAG,ICHSO,
     +      NATOMS,NGENG)
         IAORDS=IAORDS+1
         CALL CHECK(IAORDS,ISOM,ICHSO,9*NEQA,15*NEQA,JAOR(5))
       ENDIF
       IRDG=IRD(5)
       IF(IRDG.EQ.0) GO TO 8005
       DO 1665 I5=1,IRDG
          ICCX=ICCX+1
          IFL1(ICCX)=JAOR(5)
 1665  CONTINUE
 8005  CONTINUE
C
 405   CONTINUE
 300  CONTINUE
C
      IF(NGENG.EQ.0) GO TO 314
      DO 253 I=1,NGENG
         LKK1=0
         LOGO(I)=.FALSE.
         DO 251 JJ=1,ISST
            NQQ=NEQ(JJ)
            DO 251 KK=1,NQQ
               LKK1=LKK1+1
               IF(IMAGE(LKK1,I).NE.KK) GOTO 252
 251        CONTINUE
            LOGO(I)=.TRUE.
            GOTO 253
 252        NGEN=NGEN+1
 253     CONTINUE
         DO 313 IGEN=1,NGENG
            IF(LOGO(IGEN)) GOTO 313
            IGG=IGEN+1
            IF(IGEN.EQ.NGENG) GOTO 314
            DO 312 JGEN=IGG,NGENG
               IF(LOGO(JGEN)) GOTO 312
               DO 311 K=1,NATOMS
                  IF(IMAGE(K,IGEN).NE.IMAGE(K,JGEN)) GOTO 312
 311           CONTINUE
               LOGO(JGEN)=.TRUE.
               NGEN=NGEN-1
 312        CONTINUE
 313     CONTINUE
 314     CONTINUE
C
C  WRITE ARGOSIN FILE
C
      CALL WRTARG(ITAPE,TITLE,NTATOM,NSYT,NGEN,IAORDS,
     1  JJCH,IBF,NNSTV,IGR,NSTV,RR,IATM,ATNM,
     2  LOGO,IMAGE,NEQ,ICCN,NGENG,IATN)
C
C     CALCULATE THE NUMBER OF ORBITALS IN THE IRREPS
C
      LL=0
      DO 2 ITATOM=1,NTATOM
         ISHELL=NSYT(ITATOM)
         IOFFS=ICCN(ITATOM)
         DO 3 KSHELL=1,ISHELL
            NCON=SHELLSIZ(IOFFS+KSHELL)
            LL=LL+1
            IAOR=IFL1(LL)
            INF2=INFO(IAOR,2)
            IAD=INFO(IAOR,1)
            DO 4 JINF=1,INF2
               ITJ=IDUM(IAD+JINF)
               NBFS(ITJ)=NBFS(ITJ)+NCON
 4          CONTINUE
 3       CONTINUE
 2    CONTINUE
      WRITE(6,5)
 5    FORMAT(/,/,T8,'Irrep',8X,'# of functions')
      WRITE(6,6) (IBF(ISYM),NBFS(ISYM),ISYM=1,JJCH)
 6    FORMAT(T10,A3,13X,I3)
      RETURN
      END
