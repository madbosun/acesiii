
#ifdef __fortran77

C This subroutine should return the current time information
C in its parameters. Add another IFDEF group for your
C architecture/OS and enjoy!

      subroutine wallclock(year,month,day,hour,minute,second)

c#include <aces.h>

      integer year,month,day,hour,minute,second

#ifdef _AIX      
c     Please note the underscores!
      integer*4 currtime,time_
      integer*4 tarray(9)
      currtime=time_()
      call ltime_(currtime, tarray)
      second=tarray(1)
      minute=tarray(2)
      hour=tarray(3)
      day=tarray(4)
      month=tarray(5)
      year=tarray(6)+1900
c      day of week: tarray(7)
c      day of year: tarray(8)
c      daylight savings time: tarray(9)
#endif /* _AIX */

#ifdef _SUNOS
      integer*4 iarray (3)
      call idate(iarray)
      day=iarray(1)
      month=iarray(2)
      year=iarray(3)
      call itime(iarray)
      hour=iarray(1)
      minute=iarray(2)
      second=iarray(3)
#endif /* _SUNOS */

#ifdef _UNICOS
      character*8 cdate, ctime
      integer conv_char_int
      call date(cdate)
      month = 10*conv_char_int(cdate(1:1))
     &    + conv_char_int(cdate(2:2))
      day = 10*conv_char_int(cdate(4:4))
     &    + conv_char_int(cdate(5:5))
      year = 10*conv_char_int(cdate(7:7))
     &    + conv_char_int(cdate(8:8))
      call clock(ctime)
      hour = 10*conv_char_int(ctime(1:1))
     &    + conv_char_int(ctime(2:2))
      minute = 10*conv_char_int(ctime(4:4))
     &    + conv_char_int(ctime(5:5))
      second = 10*conv_char_int(ctime(7:7))
     &    + conv_char_int(ctime(8:8))
#endif /* _UNICOS */

#ifdef _IRIX
      integer*4 currtime,time_
      integer*4 tarray(9)
      currtime=time_()
      call ltime_(currtime, tarray)
      second=tarray(1)
      minute=tarray(2)
      hour=tarray(3)
      day=tarray(4)
      month=tarray(5)
      year=tarray(6)+1900
#endif /* _IRIX */

#ifdef _DEC
#   ifdef F_64BIT
      integer*4 iarray (3)
      external idate
      external itime
      call idate(iarray)
      day=iarray(1)
      month=iarray(2)
      year=iarray(3)
#   else
      integer iarray (3)
      call idate(month,day,year)
      year=year+1900
#   endif /* F_64BIT */
      call itime(iarray)
      hour=iarray(1)
      minute=iarray(2)
      second=iarray(3)
#endif /* _DEC */

c#ifdef M_IP22
c      integer iarray (3)
c      call idate (month, day, year)
c      call itime (iarray)
c      hour=iarray(1)
c      minute=iarray(2)
c      second=iarray(3)
c#endif

CC#ifdef M_HPUX
C      integer iarray (3)
C      call idate(month,day,year)
C      year=year+1900
C      call itime(iarray)
C      hour=iarray(1)
C      minute=iarray(2)
C      second=iarray(3)
CC..#ELSE
C      This subroutine needs FPP processing.
C#ENDIF

      return
      end

#elif defined(__ansi_c)

#include <time.h>

void
#ifdef C_SUFFIX
wallclock_
#else
wallclock
#endif /* C_SUFFIX */
( long * year, long * mon, long * mday, long * hour, long * min, long * sec )
{
    /*
     * Relevant members of the tm struct:
     * int tm_sec
     * int tm_min
     * int tm_hour
     * int tm_mday
     * int tm_mon
     * int tm_year
     */

    time_t    now_clock;
    struct tm now;

    time(&now_clock);
    now = *localtime(&now_clock);

 /* MUWAHAHAHAHA! ACES is not fully Y2K compliant! */
    *year = now.tm_year +1900;
    *mon  = now.tm_mon  +1;
    *mday = now.tm_mday;
    *hour = now.tm_hour;
    *min  = now.tm_min;
    *sec  = now.tm_sec;

    return;
}

#else
CRASH HARD!
#endif /* __fortran77 */
 
