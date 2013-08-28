
C THIS ROUTINE RETRIEVES A LOGICAL RECORD FROM THE JOBARC FILE.

C   INPUT:
C         LUJNK - INTEGER 
c                 > 0 missing records are reported and run aborts
c                 = 0 length of record returned in length
c                 < 0 zeros are silently returned
C         FILNAM- CHARACTER*(*) ignored
C         LABEL - IDENTIFIER FOR LOGICAL RECORD.  CHARACTER*(*) STRING.
C         LENGTH- LENGTH OF LOGICAL RECORD IN *INTEGER* WORDS.
C
C   OUTPUT:
C         ISTUFF- CONTENTS OF LOGICAL RECORD (lujnk .ne. 0 )
c         length- length of record in *integer* words ( lujnk = 0)

      SUBROUTINE GETREC(LUJNK,FILNAM,LABEL,LENGTH,ISTUFF)
      IMPLICIT INTEGER(A-Z)

      PARAMETER (LUFIL = 75)
      PARAMETER (IBUFLN = 128)
      CHARACTER*(*) LABEL
      CHARACTER*8 MARKER
      CHARACTER*80 FNAME1
      CHARACTER*(*) FILNAM
      DIMENSION ISTUFF(LENGTH),IBUF(IBUFLN)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /JOBARC/ MARKER(1000),LOC(1000),SIZE(1000),NRECS,
     &                IRECWD,IRECLN

      IF(LUJNK .GT. 0 .AND. LENGTH.EQ.0)RETURN

      CALL GFNAME('JOBARC  ',FNAME1,ILENGTH1)

          OPEN(UNIT=LUFIL,FILE=FNAME1(1:ILENGTH1),
     &         FORM='UNFORMATTED',STATUS='OLD',
     &         ACCESS='DIRECT',RECL=IRECLN)
 
c
C
C LOCATE POINTER TO THE RECORD WHICH IS TO BE READ.
C
      ILOC=LOCCHR(1000,MARKER,1,LABEL)
c
c  if record not found
c
      IF (ILOC.EQ.1001) THEN
c
         IF (LUJNK.GT.0) THEN
            WRITE(LUOUT,50)LUFIL,LABEL
 50         FORMAT(T3,'@GETREC-F, File attached to unit ',I3,
     $           ' does not contain record ',A8,'.')
            CALL ERREX
         else if ( lujnk .eq. 0 ) then
            length = -1
            close(LUFIL)
            return
         ELSE
            CALL IZERO(ISTUFF,MAX(1,LENGTH))
            close(LUFIL)
            RETURN
         ENDIF 
c
      ENDIF
c
c  record found
c
      if ( lujnk .eq. 0 ) then
         length = size(iloc)
         close(LUFIL)
         return
      endif
c     
c  this warning should not appear -- if it does not, change to ERROR
      if ( length .gt. size(iloc) ) then
         write (6,*) '@GETREC: requested length of vector ', label,
     $        ' is longer than actual'
         write (6,*) 'Actual -- ', size(iloc), ' Requested -- ', length
      endif         
C
C FIND WORD ADDRESS OF FIRST ELEMENT ON RECORD, COMPUTE ADDRESS OF LAST
C  WORD AND THEN COMPUTE BEGINNING AND ENDING RECORD NUMBERS, AND WORD
C  OFFSETS.
C
C  IWRDST - THE ABSOLUTE WORD ADDRESS WHERE THE ARRAY BEGINS.
C  IWRDEN - THE ABSOLUTE WORD ADDRESS WHERE THE ARRAY ENDS.
C  IRECST - THE PHYSICAL RECORD WHERE IWRDST IS LOCATED.
C  IRECEN - THE PHYSICAL RECORD WHERE IWRDEN IS LOCATED.
C  IOFFST - THE RELATIVE WORD ADDRESS OF IWRDST ON RECORD IRECST.
C  IOFFEN - THE RELATIVE WORD ADDRESS OF IWRDEN ON RECORD IRECEN.
C
      IWRDST=LOC(ILOC)
      IRECST=1+IWRDST/IRECWD
      IOFFST=MOD(IWRDST,IRECWD)
      IF(IOFFST.EQ.0)THEN
       IRECST=IRECST-1
       IOFFST=IRECWD
      ENDIF
      NWRDST=IRECWD-IOFFST+1
      IWRDEN=LOC(ILOC)+LENGTH-1
      IRECEN=1+IWRDEN/IRECWD
      IOFFEN=MOD(IWRDEN,IRECWD)
      IF(IOFFEN.EQ.0)THEN
       IRECEN=IRECEN-1
       IOFFEN=IRECWD
      ENDIF
      NWRDEN=IOFFEN
      IF(IRECEN.EQ.IRECST)THEN
C
C CASE ONE.  ALL INFORMATION ON ONE RECORD.
C
       READ(LUFIL,REC=IRECST,ERR=999,IOSTAT=IOS)IBUF
       CALL ICOPY(LENGTH,IBUF(IOFFST),1,ISTUFF,1)
       close(LUFIL)
       RETURN
      ELSE
C
C CASE TWO.  INFORMATION SPREAD ACROSS SEVERAL RECORDS.  FIRST GET STUFF
C  FROM FIRST RECORD.
C
       READ(LUFIL,REC=IRECST,ERR=999,IOSTAT=IOS)IBUF
       CALL ICOPY(NWRDST,IBUF(IOFFST),1,ISTUFF,1)
       IF(IRECEN.EQ.IRECST)then
         close(LUFIL)
         RETURN
       endif
C
C NOW PICK UP STUFF FROM INTERVENING RECORDS.
C
       IOFF = 1 + NWRDST
       DO 10 I=IRECST+1,IRECEN-1
        READ(LUFIL,REC=I,ERR=999,IOSTAT=IOS)IBUF
        CALL ICOPY(IRECWD,IBUF,1,ISTUFF(IOFF),1)
        IOFF=IOFF+IRECWD
10     CONTINUE
C
C NOW GET STUFF OFF LAST RECORD
C
       READ(LUFIL,REC=IRECEN,ERR=999,IOSTAT=IOS)IBUF
       CALL ICOPY(NWRDEN,IBUF,1,ISTUFF(IOFF),1)
       close(LUFIL)
       RETURN
      ENDIF
c
999   CONTINUE

      WRITE(LUOUT,*) '@GETREC: I/O error on unit ',LUFIL,
     &               ' searching for : ',LABEL
      WRITE(LUOUT,*) '         System error code : ',IOS
      CALL ERREX

      END

