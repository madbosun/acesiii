      integer max_timers
      integer max_timer_desc_len
      integer max_counter_desc_len
      integer cpu_timer
      integer elapsed_time_timer
      integer max_counters

      parameter (max_timers = 100000)
      parameter (max_counters = 100000)
      parameter (max_timer_desc_len = 40)
      parameter (max_counter_desc_len = 40)
      parameter (cpu_timer = 1)
      parameter (elapsed_time_timer = 2)

      double precision timers, tmark, timer_min, timer_max
      integer timer_type, itmark
      character*(max_timer_desc_len) tdesc
      character*(max_counter_desc_len) counter_desc
      logical do_timer

      common /timerz/timers(max_timers), 
     *               tmark(max_timers), itmark(max_timers),
     *               tdesc(max_timers), timer_type(max_timers),
     *               do_timer

      integer counters
      common /counterz/counters(max_counters), 
     *                 counter_desc(max_counters)
