
#ifdef __fortran

#   if defined(__fortran77) && defined(_UNICOS)

c PORTING:
c  - older Cray platforms (viz. SV1 generation) do not have getrusage, so
c    rename this file to timer.F

      subroutine timer(type)
      implicit none
      integer type
#include "timeinfo.com"
      double precision t
      real*8 secondr
      t = secondr()
      if (type.eq.0) then
         timein  = t
         timenow = t
         timetot = 0.0d0
         timenew = 0.0d0
      else
         timenew = t-timenow
         timenow = t
         timetot = t-timein
      end if
      return
      end

#   endif /* __fortran77 && _UNICOS */

#else /* __fortran */

#   ifdef __cplusplus
       extern "C" {
#   endif

/*
 * This initializes the f_timeinfo timer if type=0; otherwise,
 * it records the current time.
 */

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#define RUSAGE_SELF 0
#define RUSAGE_CHILDREN -1

#include "f77_name.h"
#include "f_types.h"

#include "aces.h"
#include "timeinfo.com"

void
F77_NAME(timer,TIMER)
(f_int * type)
{

    /*
     * Relevant member of the rusage struct:
     * struct timeval
     * {
     *     long tv_sec;
     *     long tv_usec;
     * } ru_utime;
     */

    struct rusage res;

    int i = getrusage(RUSAGE_SELF,&res);
    if (i == 0)
    {
        M_REAL t;
        t = res.ru_utime.tv_sec + ( res.ru_utime.tv_usec * 0.000001 );
        if (*type == 0)
        {
            f_timeinfo.timein  = t;
            f_timeinfo.timenew = 0.0;
            f_timeinfo.timetot = 0.0;
        }
        else
        {
            f_timeinfo.timenew = t - f_timeinfo.timenow;
            f_timeinfo.timetot = t - f_timeinfo.timein;
        }
        f_timeinfo.timenow = t;
    }
    else
    {
        /* error */
    }

    return;
}

#   ifdef __cplusplus
       }
#   endif

#endif /* __fortran */

