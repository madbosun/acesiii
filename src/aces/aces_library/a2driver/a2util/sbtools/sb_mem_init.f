
c This routine carves up the icore(i0) space for using setptr/relptr.
c It's clumsy, but you must communicate with this routine through
c sbcore.com by setting ineeded, dneeded, maxicore, etc.

c For dummy runs (when memknown=0), iSize sets maxmem.

c INPUT:
c int iHeap : the first integer in a heap (i.e. icore(i0) in ACES lingo)
c int iSize : the number of INTEGERS in the heap

      subroutine sb_mem_init(iHeap,iSize)
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
c      parameter (nondynimem=1000000)
c      parameter (nondyndmem=4000000)

      double precision   dcore(nondyndmem)
      common /dcore_com/ dcore
      save   /dcore_com/

      integer    icore(nondynimem)
      common / / icore

      integer             memknown,maxicore
      common /sbcore_com/ memknown,maxicore
      save   /sbcore_com/





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



      integer iHeap(*), iSize
      INTEGER c_adr, zheap, ztmp
      external c_adr

      call callstack_push('SB_MEM_INIT')

      if (iSize.le.0) then
         print *, '@SB_MEM_INIT: Assertion failed.'
         print *, '              iSize = ',iSize
         call errex
      end if

      if (memknown.ne.0) then

         if (dynmem.eq.1) then
            maxmem=iSize
            if (memknown.eq.-1) then
               ineeded=maxicore
               dneeded=(maxmem-maxicore)/iintfp
            end if
            zheap = c_adr(iHeap)

c i0 is the icore index that addresses iHeap
c shift zheap up to a 64-bit boundary
            zheap = ifltln * ((zheap+ifltln-1)/ifltln)
            ztmp = (zheap-c_adr(icore(1)))/iintln
            i0 = 1 + ztmp
            if (i0.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

c i1 is the first icore index that cannot be used (after the "icore" heap)
c (Yau: shift this up to 64-bit boundary)
            ztmp = ztmp + ineeded + iand(ineeded,iintfp-1)
            i1 = 1 + ztmp
            if (i1.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

c repeat the process for dcore
            zheap = c_adr(icore(i1))

c d0 is the dcore index that addresses zheap
            ztmp = (zheap-c_adr(dcore(1)))/ifltln
            d0 = 1 + ztmp
            if (d0.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

c d1 is the first dcore index that cannot be used (after the "dcore" heap)
            ztmp = ztmp + dneeded
            d1 = 1 + ztmp
            if (d1.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

         else
            maxmem=0
            maxicore=nondynimem
            ineeded=nondynimem
            dneeded=nondyndmem
            i0=1
            d0=1
            i1=nondynimem+1
            d1=nondyndmem+1
         end if
         iptr=i0
         dptr=d0
         call izero(icore(i0),ineeded)
         call dzero(dcore(d0),dneeded)

      else ! dummy run

         if (dynmem.eq.1) then
            maxmem=iSize
            i1=maxicore+1
            d1=(maxmem-maxicore)/iintfp+1
         else
            maxmem=0
            maxicore=nondynimem
            i1=nondynimem+1
            d1=nondyndmem+1
         end if
         i0=1
         d0=1
         iptr=1
         dptr=1
         ineeded=0
         dneeded=0

      end if

      call callstack_pop
      return
      end

