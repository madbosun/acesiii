
c This transforms a square block of a matrix
      subroutine mat_trans_sqr_sqr(in,out,inlen,outlen,inoff,outoff)


















c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler


























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      implicit none


c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/


      integer inoff,outoff,inlen,outlen
      double precision in(inlen,inlen),out(outlen,outlen)

      integer len,i

      callstack_curr='MAT_TRANS_SQR_SQR'
      len=min(inlen,outlen)
      do i=1,len
        call dcopy(len,in(inoff,inoff+i-1),1,out(outoff,outoff+i-1),1)
      enddo
      return
      end
