
c This function returns a pointer to some location in an array (either
c kscore or dcore).  It also keeps track of memory usage and memory
c requirements.
c
c The first argument (stack) is for future use.  Multiple stacks will
c be implemented for multithread processing.

      integer function setptr(stack,type,num)
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






c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end





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




c Two different "stacks" of memory will be available in the course of an
c Aces3 calculation.  Integer memory will be stored in icore.  Real memory
c will be stored in dcore.
c
c Memory allocation is done in one of two ways, dynamic or nondynamic.  This
c is controlled by the parameter dynmem.  If dynmem is 1, dynamic memory
c is used.  Otherwise, nondynamic memory is used.
c
c Two additional parameters nondynimem and nondyndmem control how much
c nondynamic memory is allocated initially in the icore and dcore stacks.
c If dynmem is 1, both of these parameters SHOULD be 1.

c Dynamic memory allocation is done in a fairly straighforward way.  The
c program runs through everything twice.  The first time is just to determine
c how much dynamic memory is required.  The second time is use to actually
c allocate the memory and use it.  In some cases, it may be difficult to
c determine how much memory to use in advance in each of the stacks.  To
c aid this, there is also the option of only running through a program a
c single time.  When this happens, both stacks are allocated initially and
c used throughout the program.  If the memory requirement for either stack
c is exceeded, the program crashes.
c
c Two parameters which are used in either case are memknown and maxicore.
c Both MUST be set in the calling program.
c
c   memknown : 0 if a dummy run is being done to determine memory usage
c              1 if the memory usage is known from a previous dummy run
c             -1 if no dummy run is done (only one run through the program)
c   maxicore : the maximum amount of icore to allocate (this MAY be
c              adjusted IF a dummy run detects that more is needed, but in
c              the case of a single run, this will be rigidly applied)
c
c This include file is required whenever calls to any of the memory functions
c (setptr and relptr) are called or when part of the execution depends on
c whether the first or second run is being performed.
c
c ***NOTE***
c It is important that icore be aligned on a floating point boundary.  One
c way which seems to insure this is to have icore the FIRST element in the
c common block.  So, make sure that no integer is ever inserted before icore
c in the common block declaration.

      integer dynmem,nondynimem,nondyndmem
      parameter (dynmem=1)
      parameter (nondynimem=1)
      parameter (nondyndmem=1)

      double precision   dcore(nondyndmem)
      common /dcore_com/ dcore
      save   /dcore_com/

      integer    kscore(nondynimem)
      common /kscore_mem / kscore

      integer            memknown,maxicore
      common /kscore_com/ memknown,maxicore
      save   /kscore_com/





c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/



      integer stack,type,num

      callstack_prev = callstack_curr
      callstack_curr = 'SETPTR'

      if (type.eq.0) then
         setptr=iptr
         iptr=iptr+num+iand(num,iintfp-1)
      else if (type.eq.1) then
         setptr=dptr
         dptr=dptr+num
      else
         write(*,*) '@SETPTR: received invalid data type'
         call errex
         stop
      end if

      if (memknown.eq.0) then
         ineeded=max(ineeded,iptr-1)
         dneeded=max(dneeded,dptr-1)
      else
         if (iptr.gt.i1 .or. dptr.gt.d1) then
            write(*,*) '@SETPTR: ran out of memory'
            call errex
            stop
         end if
      end if

      callstack_curr = callstack_prev
      return
      end

