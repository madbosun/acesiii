
c This initializes the timer if type=0; otherwise, it gets time info.

      subroutine timer(type)

#include <aces.h>
#include <const.par>
#include <timeinfo.com>

      integer type
      M_REAL t

#ifdef _UNICOS
C     Is this a system intrinsic function?
      double precision second
      t = second()
#else
      M_REAL tarray(2)
      call a3time(tarray)
      t = tarray(1)
#endif /* _UNICOS */

      if (type.eq.0) then
         timein  = t
         timenew = zero
         timetot = zero
      else
         timenew = t - timenow
         timetot = t - timein
      end if
      timenow = t

      return
      end

