
C DRIVES THE CREATION OF MOIO POINTERS FOR SYMMETRY PACKED LISTS.

      SUBROUTINE INIPCK(IRREPX,SYTYPL,SYTYPR,LIST,IARG1,IARGX2,ISET)
      IMPLICIT INTEGER (A-Z)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
#ifdef _ASSERT
      if ((SYTYPL.lt.1).or.(22.lt.SYTYPL).or.
     &    (SYTYPR.lt.1).or.(22.lt.SYTYPR)    ) then
         print *, '@INIPCK: Assertion failed.'
         print *, '   SYTYPL = ',SYTYPL
         print *, '   SYTYPR = ',SYTYPR
         print *, '   LIST   = ',LIST
         call errex
      end if
#endif /* _ASSERT */
      IARG2 = IARGX2
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,IRREPX)
         CALL UPDMOI(IRPDPD(IRREPR,SYTYPR),IRPDPD(IRREPL,SYTYPL),
     &               IRREPR,LIST,IARG1,IARG2)
         IARG1 = 0
         IARG2 = 0
      END DO
      IF (ISET.EQ.0) RETURN
      LISTX = LIST
      IF (LISTX.EQ.0) LISTX = 100
      ISYTYP(1,LISTX) = SYTYPL
      ISYTYP(2,LISTX) = SYTYPR
      RETURN
      END
