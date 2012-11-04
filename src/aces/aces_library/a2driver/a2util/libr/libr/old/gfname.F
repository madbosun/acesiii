C
C THIS SUBROUTINE GETS FOR A FILE WITH THE DEFAULT NAME
C 'NAME1' THE FULL FILENAME WHICH IS THEN RETURNED ON
C 'NAME2'. IEND2 GIVES THE LENGTH (NUMBER OF CHARACTERS)
C OF THE FILE NAME RETURNED IN 'NAME2'.
cjp
cjp for mr-bwcc moints moabcd are affixed with .iref
cjp
C
CEND
C
      SUBROUTINE GFNAME(NAME1,NAME2,IEND2)
      IMPLICIT INTEGER (A-Z)
      LOGICAL YESNO
      CHARACTER*80 NAME2,A
      CHARACTER*(*) NAME1

cjp
#include "bwcc.com"

c ----------------------------------------------------------------------

      iLen = min(8,len(Name1))
      DO 10 I = 1, iLen
         IF (NAME1(I:I).EQ.' ') THEN
            IEND = I-1
            GO TO 11
         END IF
10    CONTINUE
      IEND = iLen
11    CONTINUE

      INQUIRE(FILE='FILES',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(UNIT=90,FILE='FILES',STATUS='OLD',FORM='FORMATTED')
         REWIND(UNIT=90)
 1000    READ(90,9000,END=2000) A
 9000    FORMAT(A80)
         IF (A(1:IEND).EQ.NAME1(1:IEND).AND.
     &       A(IEND+1:IEND+1).EQ.' '        ) THEN
            DO 20 I = IEND+2, 80
               IF (A(I:I).EQ.' ') THEN
                  IEND1 = I-1
                  GO TO 21
               END IF
20          CONTINUE
            CALL ERREX
21          CONTINUE
            IEND2 = IEND1-IEND+1
            NAME2(1:IEND2) = A(IEND+2:IEND1)
            close(90)
cjp            RETURN
            goto 3000
         END IF
         GO TO 1000
 2000    CONTINUE
      END IF

      NAME2(1:IEND) = NAME1(1:IEND)
      IEND2 = IEND
      close(90)

3000  continue
      if (isbwcc) then
      if (.not.masik) then
cjp   affix with iref number for mr-bwcc
      if(name1(1:iend).eq.'MOINTS' .or. name1(1:iend).eq.'MOABCD'
     +        .or. name1(1:iend).eq.'FOCKCD') then
cjp affix number to name2
      name2(iend2+1:iend2+1)='.'
cjp ASCII assumed
      name2(iend2+2:iend2+2)=char(48+iref/10)
      name2(iend2+3:iend2+3)=char(48+iref-10*(iref/10))
      iend2=iend2+3
      endif
      endif
      endif
      if (bwgossip) write(6,*)'@GFNAME',isbwcc,masik,iref,' '
     +,name1(1:iend),'->',name2(1:iend2)

c      print *, '@GFNAME: (debug) "',NAME1,'" is "',NAME2(1:IEND2),'"'

      RETURN
      END

