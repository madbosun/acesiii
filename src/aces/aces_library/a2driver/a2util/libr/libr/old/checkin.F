
C INITIALIZES TIMING INFORMATION FOR PROGRAM

      SUBROUTINE CHECKIN
      IMPLICIT INTEGER (A-Z)

#include "bwcc.com"
#include "timeinfo.com"

cjp to make the timer consistent in MR
      if ((.not.isbwcc).or.(iref.eq.myproc)) call zclock(0)

      RETURN
      END

