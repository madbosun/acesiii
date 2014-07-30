
c This function returns the amount of memory remaining.

      integer function memleft(stack,type)
      implicit none
















c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler


























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c i0,i1    : The first element in icore that may be used and the first element
c            that may not be used.
c iptr     : The pointer to the current location in icore
c ineeded  : The size of icore needed for the calculation (determined by the
c            initial run).  This is only used when dynamic memory is used.
c
c In a dummy run, i0 starts at 1, iptr is incremented, and ineeded goes from
c 0 to ineeded.  i1 is initially set at maxicore.  These values are used if
c nondynamic memory is used.  Note that if dynamic memory is used, ineeded
c MAY exceed i1.  In this case (provided there dcore has not used all the
c remaining memory), the amount of memory given to icore is increased.  On
c the other hand, if icore does not need all that it was given, it gives
c the remaining to dcore.
c
c When memory is actually allocated, i0 is set to the beginning element,
c i1 is set to i0+ineeded, and iptr ranges from i0 to i1.
c
c d0,dpt,dneeded,d1
c          : Similar for dcore
c
c maxmem   : Maximum memory to allocate (set in ZMAT file)

      integer             i0,iptr,ineeded,i1,d0,dptr,dneeded,d1,maxmem
      common /ks_mem_com/ i0,iptr,ineeded,i1,d0,dptr,dneeded,d1,maxmem
      save   /ks_mem_com/




c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/



      integer stack,type

      callstack_prev = callstack_curr
      callstack_curr='MEMLEFT'

      if (type.eq.0) then
         memleft=i1-iptr
      else if (type.eq.1) then
         memleft=d1-dptr
      else
         write(*,*) '@MEMLEFT: received invalid data type'
         call errex
         stop
      end if

      callstack_curr = callstack_prev
      return
      end

