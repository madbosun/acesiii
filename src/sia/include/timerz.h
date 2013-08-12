      integer max_timers
      integer max_allocated_timers
      integer current_last_timer
      integer max_timer_desc_len
      integer cpu_timer
      integer elapsed_time_timer
      integer times_exec_timer
      integer average_unit_timer

      parameter (max_allocated_timers = 100000)
      parameter (max_timer_desc_len = 40)
      parameter (cpu_timer = 1)
      parameter (elapsed_time_timer = 2)
      parameter (times_exec_timer = 3)
      parameter (average_unit_timer = 4)

      double precision timers, tmark, timer_min, timer_max
      double precision timer_optl, timer_ovrhead, interim_timers
      double precision interim_start_timers
      integer timer_type, itmark, timer_times
      character*(max_timer_desc_len) tdesc
      logical do_timer

      common /timerz/tdesc(max_allocated_timers),
     *               timers(max_allocated_timers),
     *               interim_timers(max_allocated_timers),
     *               interim_start_timers(max_allocated_timers),
     *               tmark(max_allocated_timers),
     *               itmark(max_allocated_timers),
     *               timer_type(max_allocated_timers),
     *               timer_times(max_allocated_timers),
     *               timer_optl, timer_ovrhead,
     *               current_last_timer, do_timer,
     *               max_timers
