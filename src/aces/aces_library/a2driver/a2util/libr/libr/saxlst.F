      SUBROUTINE SAXLST(ICORE,MAXCOR,LSTINA,LSTINB,LSTOUT,FA,FB)
C
C THIS ROUTINE READS TWO LISTS, A AND B, AND FORMS A THIRD
C (C) WHICH IS EQUAL TO FA * A + FB * B.  LISTS A AND B ARE
C READ FROM LSTINA AND LSTINB AND THE RESULTANT C IS WRITTEN
C TO LSTOUT.  THE SIZES ARE ASSUMED TO BE THE SAME (NSIZE).
C
C A couple of silly bugs corrected by RPM 5/13/94.
C
CEND
      IMPLICIT INTEGER (A-H,O-Z)
      DOUBLE PRECISION ONE,FA,FB
      character *(*)rtnn,rtn
      PARAMETER(ONE=1.0D+00)
      parameter (rtnn='saxlst', rtn='  @'//rtnn)
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
C DETERMINE SIZES AND SEE IF THIS WILL GO WITH PUTALL AND GETALL.
C
      LSTTYP=LSTINA
      IF(LSTINA.GT.100)LSTTYP=LSTTYP
      NSIZE=ISYMSZ(ISYTYP(1,LSTTYP),ISYTYP(2,LSTTYP))
C      IF(NSIZE.LT.MAXCOR)THEN !!!!!!  No no no RPM
      IF(2*NSIZE .LT. MAXCOR/IINTFP)THEN
C         write (*,*) rtn,'-I: getall in core'
         I000=1
         I010=I000+IINTFP*NSIZE
         I020=I010+IINTFP*NSIZE
         if ( i020 .gt. maxcor ) then
            write (*,*) rtn,'-F: no room for both lists'
            call insmem(rtnn, i020, maxcor)
         endif
         CALL GETALL(ICORE(I000),NSIZE,1,LSTINA)
         CALL SSCAL(NSIZE,FA,ICORE(I000),1)
         CALL GETALL(ICORE(I010),NSIZE,1,LSTINB)
         CALL SSCAL(NSIZE,FB,ICORE(I010),1)
         CALL SAXPY(NSIZE,ONE,ICORE(I000),1,ICORE(I010),1)
         CALL PUTALL(ICORE(I010),NSIZE,1,LSTOUT)
      ELSE
C         write (*,*) 'out of core'
C         write (*,*)'  pair type of list ', lsttyp,
C     &     ' is (',  isytyp(1,lsttyp),',', isytyp(2,lsttyp),')'
         DO 10 IRREP=1,NIRREP
C            write (*,*)rtn,'-I: irrep ', irrep
            DISSIZ=IRPDPD(IRREP,ISYTYP(1,LSTTYP))
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LSTTYP))
C
            NSIZE=DISSIZ*NUMDIS
            I000=1
            I010=I000+IINTFP*NSIZE
            I020=I010+IINTFP*NSIZE
            IF(I020.GT.MAXCOR) then
               write (*,*) rtn,'-F: no room for out of core ',
     &           'irrep ', irrep
               CALL  INSMEM('SAXLST',I020,MAXCOR)
            endif
            CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREP,LSTINA)
            CALL SSCAL(NSIZE,FA,ICORE(I000),1)
            CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LSTINB)
            CALL SSCAL(NSIZE,FB,ICORE(I010),1)
            CALL SAXPY(NSIZE,ONE,ICORE(I000),1,ICORE(I010),1)
            CALL PUTLST(ICORE(I010),1,NUMDIS,1,IRREP,LSTOUT)
 10      CONTINUE
      ENDIF
C
      RETURN
      END
