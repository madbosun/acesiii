      SUBROUTINE SPINAD3(IRREP,NUM1,DIS,NUM,A,SCR,ISCR)
C
C THIS ROUTINE SPIN ADAPTS A MATRIX A BY
C
C  A(J,I,K,L) ----> 2*A(J,I,K,L)-A(I,J,K,L) 
C
C NOTE THAT THE SPIN ADAPTION IS CARRIED OUT HERE IN PLACE AND
C ONLY TWO ADDITIONAL ARRAYS OF SIZE NUM ARE REQUIRED
C
C INPUT : IRREP .......... IRREP OF THE GIVEN PAR OF A
C         NUM1 ........... POPULATION VECTOR OF I AND J
C         DIS ............ DISTRIBUTION SIZE OF A
C         NUM ............ NUMBER OF DISTRIBUTIONS IN A
C         A .............. INPUT MATRIX A
C         SCR,ISCR ....... TWO SCRATCH ARRAYS OF DIMENSION NUM
C
C  OUTPUT : A ............ SPIN ADAPTED MATRIX A
C 
CEND
C
C CODED SEPTEMBER/90  JG
C
C The SCR array is no longer needed, but the reference is kept for
C   compatibility.
C SG 6/95
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS,DIRPRD
      DIMENSION A(DIS,NUM),SCR(1),ISCR(2),NUM1(8),IP(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
C
      DATA TWO /2.0D+00/
C
      IP(1)=0
      DO 10 IRREPJ=1,NIRREP-1
       IRREPI=DIRPRD(IRREP,IRREPJ)
       IP(IRREPJ+1)=IP(IRREPJ)+NUM1(IRREPJ)*NUM1(IRREPI)
10    CONTINUE
       
C
C  GET FIRST THE NEW ADDRESSES AND STORE THEM IN ISCR
C
C  The addresses will be stored in ISCR in the order IJ(1), IJTR(1),
C    IJ(2), IJTR(2), ...  So, for each IJ, its transpose will be the next
C    location in ISCR.  Each memory location will appear once at most.
C
      INDX = 0
      DO 20  IRREPJ=1,NIRREP
        NUMJ=NUM1(IRREPJ)
        IRREPI=DIRPRD(IRREP,IRREPJ)
        NUMI=NUM1(IRREPI)
        DO 16 J=1,NUMJ
          DO 15 I=1,NUMI
            IND1=IP(IRREPJ)+(J-1)*NUMI+I
            IND2=IP(IRREPI)+(I-1)*NUMJ+J
            IF (IND1 .LT. IND2) THEN
              ISCR(INDX+1) = IND1
              ISCR(INDX+2) = IND2
              INDX = INDX + 2
            ENDIF
 15       CONTINUE
 16     CONTINUE
 20   CONTINUE
C
C  NOW SPIN ADAPT A
C  SPIN ADAPT A(IJ,L) and A(IJTR,L)
C
      DO 100 L = 1,NUM
CDIR$ IVDEP
*VOCL LOOP,NOVREC
        DO 5 M = 1,INDX,2
          IJ = ISCR(M)
          IJTR = ISCR(M+1)
          TEMP = A(IJ,L)
          A(IJ,L) = TWO*A(IJ,L) - A(IJTR,L)
          A(IJTR,L) = TWO*A(IJTR,L) - TEMP
 5      CONTINUE
C
C ALL DONE FOR IJ AND IJTR
C
 100  CONTINUE 
C
      RETURN
      END
