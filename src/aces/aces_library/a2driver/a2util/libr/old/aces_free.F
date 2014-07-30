
      SUBROUTINE ACES_FREE(FREWRD,IFIRST,IPOS)
      IMPLICIT INTEGER (A-Z)

      DOUBLE PRECISION TELAP

#include "machsp.com"

#ifdef _UNICOS
      IFRE=IPOS+LOC(IFIRST)-1
      ILEN=IHPLEN(IFRE,IERR,0)
      CALL HPDEALLC(IFRE,IERR,0)
#else
      IFRE=lshift(IPOS,2)+LOC(IFIRST)-4
      CALL C_FREE(IFRE)
#endif /* _UNICOS */

      write(*,*) '@ACES_FREE: Released ',FREWRD,
     &           ' integer words of core memory.'

      RETURN
      END

