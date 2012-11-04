
C DETERMINES TOTAL EXECUTION TIME FOR PROGRAM

      SUBROUTINE CHECKOUT
      IMPLICIT INTEGER (A-Z)

#include "bwcc.com"
#include "timeinfo.com"

cjp - only after we are finished with the last reference in MR
      if ((.not.isbwcc).or.(iref+nproc.gt.nref)) call zclock(1)

      RETURN
      END

