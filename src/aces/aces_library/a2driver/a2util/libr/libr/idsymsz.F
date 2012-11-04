
C THIS ROUTINE RETURNS THE TOTAL SIZE OF A SYMMETRY PACKED
C DISTRIBUTION (FOR ALL IRREPS).  ITYPL IS THE SYMMETRY TYPE
C OF THE DISTRIBUTION MEMBERS, ITYPR IS THE SYMMETRY TYPE OF
C THE DISTRIBUTIONS.  NOTE THAT CHANGING THE MEANING OF THESE
C VARIABLES MAKES NO DIFFERENCE.

      INTEGER FUNCTION IDSYMSZ(IRREPX,ITYPL,ITYPR)
      IMPLICIT INTEGER (A-Z)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                 DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
#ifdef _ASSERT
      if ((IRREPX.lt.1).or.(NIRREP.lt.IRREPX)) then
         print *, '@IDSYMSZ: Assertion failed.'
         print *, '   IRREPX = ',IRREPX
         print *, '   NIRREP = ',NIRREP
         call errex
      end if
      if ((ITYPL.lt.1).or.(22.lt.ITYPL).or.
     &    (ITYPR.lt.1).or.(22.lt.ITYPR)    ) then
         print *, '@IDSYMSZ: Assertion failed.'
         print *, '   ITYPL = ',ITYPL
         print *, '   ITYPR = ',ITYPR
         call errex
      end if
#endif /* _ASSERT */ 
      IDSYMSZ = 0
      DO IRREPR = 1, NIRREP 
         IRREPL = DIRPRD(IRREPR,IRREPX)
         IDSYMSZ =   IDSYMSZ
     &             + ( IRPDPD(IRREPL,ITYPL)*IRPDPD(IRREPR,ITYPR) )
      END DO
      RETURN
      END
