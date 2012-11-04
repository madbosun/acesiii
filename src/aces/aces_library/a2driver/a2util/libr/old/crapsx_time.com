#ifndef _CRAPSX_TIME_COM_
#define _CRAPSX_TIME_COM_
c crapsx_time.com : begin
      integer year_i, mon_i, mday_i, hour_i, min_i, sec_i
      integer year_o, mon_o, mday_o, hour_o, min_o, sec_o
      common /crapsx_time/
     &        year_i, mon_i, mday_i, hour_i, min_i, sec_i,
     &        year_o, mon_o, mday_o, hour_o, min_o, sec_o
c crapsx_time.com : end
#endif /* _CRAPSX_TIME_COM_ */
