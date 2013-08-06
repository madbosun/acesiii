
c This routine switches off the timing statistics that are printed by
c aces_fin.

      subroutine aces_time_hush
      implicit none
c aces_time.com : begin

c These six times hold the timing data from aces_init (*_in) to
c aces_fin (*_out). utime and stime count the total number of user
c and system seconds since the start of the process. rtime counts
c the total number of real seconds since 1 January 1970. ame_timed
c is a logical flag set in aces_init telling aces_fin to print out
c a timing summary.

      external aces_bd_aces_time

      double precision   ame_utime_in,  ame_stime_in,  ame_rtime_in,
     &                   ame_utime_out, ame_stime_out, ame_rtime_out
      logical            ame_timed
      common /aces_time/ ame_utime_in,  ame_stime_in,  ame_rtime_in,
     &                   ame_utime_out, ame_stime_out, ame_rtime_out,
     &                   ame_timed
      save   /aces_time/

c aces_time.com : end
      ame_timed = .false.
      return
      end

