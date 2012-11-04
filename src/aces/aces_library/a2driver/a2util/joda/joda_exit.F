      subroutine joda_exit(istat, string)
      integer istat
      character*(*) string

      if (istat.ne.0) print *,string,': status = ',istat
      stop
      end
