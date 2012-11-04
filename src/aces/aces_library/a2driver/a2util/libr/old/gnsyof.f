      SUBROUTINE GNSYOF
C
C THIS ROUTINE CREATES SYMMETRY OFFSET VECTORS
C
C THE NUMBERING SCHEME IS AS FOLLOWS :
C
C  a<b    (alpha) [1]
C  a<b    (beta)  [2]
C  i<j    (alpha) [3]
C  i<j    (beta)  [4]
C  a<=b   (alpha) [5]
C  a<=b   (beta)  [6]
C  i<=j   (alpha) [7]
C  i<=j   (beta)  [8]
C  a,i    (alpha) [9]
C  a,i    (beta)  [10]
C  a,i    (AB)    [11]
C  a,i    (BA)    [12]
C  a,b    (AB)    [13]
C  i,j    (AB)    [14]
C  a,b    (AB)    [15]
C  i,a    (alpha) [16]
C  i,a    (beta)  [17]
C  i,a    (AB)    [18]
C  a,b    (alpha) [19]
C  a,b    (beta)  [20]
C  i,j    (alpha) [21]
C  i,j    (beta)  [22]
C  a,b    (BA)    [23]
C  i,j    (BA)    [24]
C  i,a    (BA)    [25]
C 
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*3 TYPE(2,25)
      CHARACTER*4 PACK(25)
      DIMENSION ISYM(8,8),SPIN(2,25),POPL(8),POPR(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA TYPE /'VRT','VRT' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','OCC' , 'VRT','VRT' , 'VRT','VRT' ,
     &           'OCC','OCC' , 'OCC','OCC' , 'VRT','OCC' ,
     &           'VRT','OCC' , 'VRT','OCC' , 'VRT','OCC' ,
     &           'VRT','VRT' , 'OCC','OCC' , 'VRT','VRT' ,
     &           'OCC','VRT' , 'OCC','VRT' , 'OCC','VRT' ,
     &           'VRT','VRT' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','OCC' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','VRT' /
      DATA SPIN / 1,1 , 2,2 , 1,1 , 2,2 , 1,1 , 2,2 , 1,1 , 2,2 ,
     &            1,1 , 2,2 ,
     &            1,2 , 2,1 , 1,2 , 1,2 , 1,2 , 1,1 , 2,2 , 1,2 ,
     &            1,1 , 2,2 ,
     &            1,1 , 2,2 , 2,1 , 2,1 , 2,1 / 
      DATA PACK /'PACK','PACK','PACK','PACK','PCK2',
     &           'PCK2','PCK2','PCK2','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL'/
C
      NNM1O2(I)=(I*(I-1))/2
      NNP1O2(I)=(I*(I+1))/2
C
C CHECK TO MAKE SURE POPULATION COUNTERS HAVE BEEN INITIALIZED
C
      NTOT=0
      DO 1 ISPIN=1,2
       DO 2 IRREP=1,8
        NTOT=NTOT+POP(IRREP,ISPIN)
2      CONTINUE
1     CONTINUE 
C
      IF(NTOT.EQ.0)RETURN
C
      DO 5 ITYPE=1,25
       IF(TYPE(1,ITYPE).EQ.'VRT'.AND.SPIN(1,ITYPE).EQ.1)THEN
        CALL ICOPY(8,VRT(1,1),1,POPL,1)
       ELSEIF(TYPE(1,ITYPE).EQ.'VRT'.AND.SPIN(1,ITYPE).EQ.2)THEN
        CALL ICOPY(8,VRT(1,2),1,POPL,1)
       ELSEIF(TYPE(1,ITYPE).EQ.'OCC'.AND.SPIN(1,ITYPE).EQ.1)THEN
        CALL ICOPY(8,POP(1,1),1,POPL,1)
       ELSEIF(TYPE(1,ITYPE).EQ.'OCC'.AND.SPIN(1,ITYPE).EQ.2)THEN
        CALL ICOPY(8,POP(1,2),1,POPL,1)
       ENDIF
       IF(TYPE(2,ITYPE).EQ.'VRT'.AND.SPIN(2,ITYPE).EQ.1)THEN
        CALL ICOPY(8,VRT(1,1),1,POPR,1)
       ELSEIF(TYPE(2,ITYPE).EQ.'VRT'.AND.SPIN(2,ITYPE).EQ.2)THEN
        CALL ICOPY(8,VRT(1,2),1,POPR,1)
       ELSEIF(TYPE(2,ITYPE).EQ.'OCC'.AND.SPIN(2,ITYPE).EQ.1)THEN
        CALL ICOPY(8,POP(1,1),1,POPR,1)
       ELSEIF(TYPE(2,ITYPE).EQ.'OCC'.AND.SPIN(2,ITYPE).EQ.2)THEN
        CALL ICOPY(8,POP(1,2),1,POPR,1)
       ENDIF
C        
       DO 10 IRRDPD=1,NIRREP
        ILOC=1
C
C FULL STORAGE MODE
C
        IF(PACK(ITYPE).EQ.'FULL')THEN
         DO 20 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRRDPD)
          ISYMOFF(IRREPR,IRRDPD,ITYPE)=ILOC
          ILOC=ILOC+POPL(IRREPL)*POPR(IRREPR)
20       CONTINUE
        ELSEIF(PACK(ITYPE).EQ.'PACK')THEN
         DO 30 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRRDPD)
          IF(IRREPL.LT.IRREPR)THEN
           ISYMOFF(IRREPR,IRRDPD,ITYPE)=ILOC
           ILOC=ILOC+POPL(IRREPL)*POPR(IRREPR)
          ELSEIF(IRREPL.EQ.IRREPR)THEN
           ISYMOFF(IRREPR,IRRDPD,ITYPE)=ILOC
           ILOC=ILOC+NNM1O2(POPL(IRREPL))
          ENDIF
30       CONTINUE
        ELSEIF(PACK(ITYPE).EQ.'PCK2')THEN
         DO 40 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRRDPD)
          IF(IRREPL.LT.IRREPR)THEN
           ISYMOFF(IRREPR,IRRDPD,ITYPE)=ILOC
           ILOC=ILOC+POPL(IRREPL)*POPR(IRREPR)
          ELSEIF(IRREPL.EQ.IRREPR)THEN
           ISYMOFF(IRREPR,IRRDPD,ITYPE)=ILOC
           ILOC=ILOC+NNP1O2(POPL(IRREPL))
          ENDIF
40       CONTINUE
        ENDIF
C
10     CONTINUE
5     CONTINUE
C
      RETURN
      END
