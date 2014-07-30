
c This routine prints the elapsed time.

c INPUT
c    integer i : behavior flag
c              = 0 ; record the wall clock for future use
c             != 0 ; print the elapsed time since zclock(0) was called

#define PRINT_ZCLOCK_INFO
c#define _DEBUG_ZCLOCK

      subroutine zclock(i)
      implicit none

#include "aces.h"
#include "timeinfo.com"

      integer i

      double precision sc,                 percent
      integer          hc, mc, hw, mw, sw, wall
      integer itmp

      integer days(12), day0, day
      data days /0,31,59,90,120,151,181,212,243,273,304,334/

      logical set_up
      data set_up /.false./

      integer             y0,mo0,d0,h0,mi0,s0
      double precision                        cpu0
      common /zclock_com/ y0,mo0,d0,h0,mi0,s0,cpu0
      save   /zclock_com/

      integer y, mo, d, h, mi, s
      double precision cpu

c ----------------------------------------------------------------------

      if (i.eq.0) then
         call wallclock(y0,mo0,d0,h0,mi0,s0)
         call timer(1)
         cpu0 = timenow
         set_up = .true.
         return
      end if

      if (.not.set_up) then
         write(*,*) '@ZCLOCK: zclock called before proper set up'
         return
      end if

      call wallclock(y,mo,d,h,mi,s)
      call timer(1)

#ifdef _DEBUG_ZCLOCK
      write(*,*) '@ZCLOCK: cpu0 = ',cpu0,' and timenow = ',timenow
      write(*,*) '         wall0 = ',y0,mo0,d0,h0,mi0,s0
      write(*,*) '         wall1 = ',y, mo, d, h, mi, s
#endif

c   o SYSTEM TIME
      cpu  = timenow - cpu0
      hc   = cpu / 3600.0
      itmp = -3600 * hc
      mc   = ( cpu + itmp ) / 60.0
      sc   = cpu + itmp + ( -60 * mc )

c   o USER TIME
      if ((mo0.gt.2).and.(mod(y0,4).eq.0)) then
         day0 = days(mo0) + d0 + 1
      else
         day0 = days(mo0) + d0
      end if
      if ((mo .gt.2).and.(mod(y ,4).eq.0)) then
         day = days(mo ) + d  + 1
      else
         day = days(mo ) + d
      end if
c     Calculations longer than one year are improbable and not supported.
      if (y0.ne.y) day = day + 365

      wall =   (day - day0 )*86400
     &       + (h   - h0   )*3600
     &       + (mi  - mi0  )*60
     &       + (s   - s0   )
      hw   = wall / 3600
      itmp = -3600 * hw
      mw   = ( wall + itmp ) / 60
      sw   = wall + itmp + ( -60 * mw )

      if (wall.gt.0) then
         percent = ( cpu / wall ) * 100.0
      else
         percent = 100.0
      endif

#ifdef PRINT_ZCLOCK_INFO
c      print 9000,hc,mc,sc,hw,mw,sw,percent
c 9000    format('  cpu time:', I4,' h ',I2,' min ',F5.2,' sec;  ',
c     &   'real time:', I4,' h ',I2,' min ',I2,' sec  (',F5.1,'%)'/)
c        CPU TIME: 999 H 23 MIN 65.32 SEC   REAL TIME: 999 H 23 MIN 76
c        SEC  (100.0%)

      itmp = cpu
      write(*,9000) '  AME TIMING:  ',
     &   'cpu ', hc,':',mc,':',sc,' (',itmp,') ;  ',
     &   'real ',hw,':',mw,':',sw,' (',wall,')'
 9000 format(a,
     &       1(a,i4,a,i2,a,f5.2,a,i8,a),
     &       1(a,i4,a,i2,a,i2  ,a,i8,a)
     &      )

#endif /* PRINT_ZCLOCK_INFO */

      return
      end

