
/*
 * FORTRAN callable routine to get the number of seconds the process has been
 * running
 */

#include <sys/times.h>
/* #include "aces.h" */
#define M_REAL double

struct zts
{
    M_REAL usr;
    M_REAL sys;
};

void
#ifdef _UNICOS
     A3TIME
#else
#  ifdef C_SUFFIX
     a3time_
#  else
     a3time
#  endif /* C_SUFFIX */
#endif /* _UNICOS */
(struct zts * secs)
{
    struct tms buf;

    times(&buf);

    /*
     * tms_utime := CPU time charged for the execution of user instructions.
     * tms_stime := CPU time charged for execution by the system on behalf of
     *              the process.
     */

    secs->usr = buf.tms_utime / 100.0;
    secs->sys = buf.tms_stime / 100.0;
}

