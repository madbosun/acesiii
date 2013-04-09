C
C DRIVES THE CREATION OF MOIO POINTERS FOR SYMMETRY PACKED LISTS.
C
      SUBROUTINE INIPCK3(SYTYPL,SYTYPR,LIST,IARG1,IARGX2)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*6 SYMLST(22)
      DIMENSION POPLFT(8),POPRHT(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      DATA SYMLST /'SVAVA0','SVBVB0','SOAOA0','SOBOB0',
     &             'SVAVA1','SVBVB1','SOAOA1','SOBOB1',
     &             'SVAOA2','SVBOB2','SOBVA2','SVBOA2',
     &             'SVAVB2','SOAOB2','SVAVB2','SOAVA2',
     &             'SOBVB2','SOAVB2','SVAVA2','SVBVB2',
     &             'SOAOA2','SOBOB2'/
#ifdef _ASSERT
      if ((SYTYPL.lt.1).or.(22.lt.SYTYPL).or.
     &    (SYTYPR.lt.1).or.(22.lt.SYTYPR)    ) then
         print *, '@INIPCK3: Assertion failed.'
         print *, '   SYTYPL = ',SYTYPL
         print *, '   SYTYPR = ',SYTYPR
         print *, '   LIST   = ',LIST
         call errex
      end if
#endif /* _ASSERT */ 
      IARG2=IARGX2
      CALL GETREC(20,'JOBARC',SYMLST(SYTYPL)//'X ',NIRREP,POPLFT)
      CALL GETREC(20,'JOBARC',SYMLST(SYTYPR)//'X ',NIRREP,POPRHT)
      DO 10 IRREP=1,NIRREP
       CALL UPDMOI(POPRHT(IRREP),POPLFT(IRREP),IRREP,LIST,IARG1,IARG2)
       IARG1=0
       IARG2=0
C
C      This is not a good thing to do now that ISYTYP handles up to
C      500 lists!     db 920531
C      IF(LIST.GT.100) GO TO 10
C
       listx=list
       if(listx.eq.0)listx=100
       ISYTYP(1,LISTx)=SYTYPL
       ISYTYP(2,LISTx)=SYTYPR
10    CONTINUE
      RETURN
      END
