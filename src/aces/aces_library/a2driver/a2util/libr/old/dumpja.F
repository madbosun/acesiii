
C THIS ROUTINE SERVES AS AN INTERFACE BETWEEN THE GETREC/PUTREC I/O
C TO THE JOBARC FILE AND THE DIFFERENT PROGRAMS IN THE PROGRAM SYSTEM.
C IT IS A DEPENDENT OF BOTH CRAPSI AND CRAPSO AND SHOULD NEVER NEED
C TO BE CALLED EXPLICITLY.
c
c Now opens or initializes JOBARC, so GET/PUTREC do not have to.
C
C INPUT:
C        IENTRY - 'I' FOR INITIATION OF PROGRAM (READS IN WORD ADDRESSES
C                         AND LABELS FROM JASMRY)
C
C                 'O' FOR EXIT FROM PROGRAM (DUMPS WORD ADDRESS AND LABELS
C                         TO JASMRY)

      SUBROUTINE DUMPJA(IENTRY)
      IMPLICIT INTEGER (A-Z)
      LOGICAL YESNO, isopn
      CHARACTER*1 IENTRY
      CHARACTER*8 MARKER
      CHARACTER*80 FNAME1,FNAME2
      integer lufil
      parameter ( lufil = 75 )
      PARAMETER (IBUFLN = 128)
      integer ibuf(ibufln)
      PARAMETER (LENGTH = 1000)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
 
#include "icdacc.com"

      COMMON /JOBARC/ MARKER(LENGTH),LOC(LENGTH),SIZE(LENGTH),NRECS,
     &                IRECWD,IRECLN

      CALL GFNAME('JOBARC  ',FNAME1,ILENGTH1)
      CALL GFNAME('JAINDX  ',FNAME2,ILENGTH2)

      IF(IENTRY.EQ.'O')THEN
c         close ( lufil )
         OPEN(UNIT=lufil,FILE=FNAME2(1:ILENGTH2),FORM='UNFORMATTED',
     &        STATUS='UNKNOWN')
         REWIND(lufil)
         WRITE(lufil)MARKER,LOC,SIZE,NRECS
         CLOSE(UNIT=lufil,STATUS='KEEP')
      ELSE
         INQUIRE(FILE=FNAME2(1:ILENGTH2),EXIST=YESNO)
         IF(YESNO)THEN
            OPEN(UNIT=lufil,FILE=FNAME2(1:ILENGTH2),
     &           FORM='UNFORMATTED',STATUS='OLD')
            REWIND(lufil)
            READ(lufil)MARKER,LOC,SIZE,NRECS
            CLOSE(UNIT=lufil,STATUS='KEEP')
         ENDIF
         INQUIRE(FILE=FNAME1(1:ILENGTH1),EXIST=YESNO,OPENED=ISOPN)
         IF(.NOT.YESNO)THEN
            IRECWD=IBUFLN
c YAU - old
cC Nevin
c            if(iintfp.eq.1)then
c               ITMPP=IINTLN
c            else
c               ITMPP=idaccm
c            endif
c            IRECLN=IBUFLN*ITMPP
c YAU - new
            IRECLN=IBUFLN*idaccm
c YAU - end
            nrecs = 0
C Bug fix, AP 03/25/97.
            CALL IZERO(IBUF,IBUFLN)
            OPEN(UNIT=LUFIL,FILE=FNAME1(1:ILENGTH1),
     &           FORM='UNFORMATTED',STATUS='NEW',
     &           ACCESS='DIRECT',RECL=IRECLN)
            WRITE(LUFIL,REC=1)IBUF
            close(lufil)
            DO 5 I=1,1000
               MARKER(I)='OPENSLOT'
 5          CONTINUE
         ENDIF
         IF(.NOT.ISOPN)THEN
            IRECWD=IBUFLN
c YAU - old
cC Nevin
c            if(iintfp.eq.1)then
c               ITMPP=IINTLN
c            else
c               ITMPP=idaccm
c            endif
c            IRECLN=IBUFLN*ITMPP
c YAU - new
            IRECLN=IBUFLN*idaccm
c YAU - end
c            OPEN(UNIT=LUFIL,FILE=FNAME1(1:ILENGTH1),
c     &           FORM='UNFORMATTED',STATUS='OLD',
c     &           ACCESS='DIRECT',RECL=IRECLN)
         ENDIF
      ENDIF

      RETURN
      END
