
c This routine initializes the ACES environment for use by any program.
c It or aces_init_rte should be called before any other ACES routine.

c INPUT
c logical BALLOCMEM : do (not) allocate memory according to IFLAGS(h_IFLAGS_mem)
c                     NOTE: If memory is allocated, then aces_init will
c                           initialize the I/O subsystem and the cache
c                           subsystem.

c INPUT/OUTPUT
c int ICORE(*) : a reference integer (address) used to anchor a malloc'ed heap,
c                the value of which is set to the byte address of the heap
c                (not valid on 32-bit builds with 64-bit pointers)

c OUTPUT
c int ICORENDX : the integer index that addresses the malloc'ed heap
c int ICOREDIM : the number of addressable integers at ICORE(ICORENDX)
c int IUHF : the spin state of the reference wavefunction
c            = 0; spin-  restricted (RHF)
c            = 1; spin-unrestricted (UHF)
















c#define _DEBUG_ACES_INIT

      subroutine aces_init(iCore,iCoreNdx,iCoreDim,iUHF,bAllocMem)
      implicit none

c ARGUMENTS
      integer iCore(*), iCoreNdx, iCoreDim, iUHF
      logical bAllocMem

c EXTERNAL FUNCTIONS
      INTEGER c_adr
      external c_adr

c PARAMETERS
      integer iMemMin, iMemInc
      parameter (iMemMin=1*1024*1024,iMemInc=1*1024*1024)

c INTERNAL VARIABLES
      double precision dTmp
      integer i0, iMem, iTmp

c TIMING VARIABLES
      integer stime_sec, stime_usec
      integer utime_sec, utime_usec
      integer rtime_sec, rtime_usec

c COMMON BLOCKS
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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
c BWCC needs aces_init to be reentrant, so an initialized common value must
c exist. The blockdata routine is appended to the end of this file.
      external aces_bd_aces_reflag
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/

c ----------------------------------------------------------------------


c   o initialize timing information
      call c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      call c_gtod   (rtime_sec,rtime_usec)
      ame_stime_in = (1.d-6 * stime_usec) + stime_sec
      ame_utime_in = (1.d-6 * utime_usec) + utime_sec
      ame_rtime_in = (1.d-6 * rtime_usec) + rtime_sec
      ame_timed = .true.


c   o setup stdout buffering to see output immediately
      call bufferlines

c   o initialize the runtime environment (held in common blocks)
      call aces_init_rte

c   o gather parallel statistics (needed by gfname in ja_init)
      call aces_com_parallel_aces

c   o initialize the job archive subsystem (getrec/putrec)
      call aces_ja_init
c???      call getrec(1,'JOBARC','IENDSTAT',1,iStat)

c   o load the ACES State Variables
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)

c   o allocate memory and initialize the I/O subsystem
      if (bAllocMem) then

c      o load the requested core size
         iMem = iflags(36)

c      o allocate core memory
         if (xFlag.eq.0) then
            iCore(1) = 0
            do while ((iCore(1).eq.0).and.(iMem.gt.iMemMin))
               call aces_malloc(iMem,iCore,i0)
               if (iCore(1).eq.0) iMem = iMem - iMemInc
            end do
            if (iMem.lt.iflags(36)) then
               print *, '@ACES_INIT: MEMORY WARNING!'
               print *, '            requested ',iflags(36),' integers'
               print *, '            allocated ',iMem,' integers'
            end if
            if (iCore(1).ne.0) then
               xFlag = 1
            else
c               print *, '@ACES_INIT: Request for ',iMem,
c     &                  ' integers of memory failed.'
               print *, '@ACES_INIT: unable to allocate at least ',
     &                  iMemMin,' integers of memory'
               call aces_exit(1)
            end if
         end if

c      o initialize the I/O subsystem ('T' creates and initializes the cache)
         call aces_io_init(iCore,i0,iMem,.true.)

c      o transfer the iCore statistics
         iCoreNdx = i0
         iCoreDim = iMem

c     else if (.not.bAllocMem)
      else

         iCore(1) = 0
         iCoreNdx = 1
         iCoreDim = 1

c     end if (bAllocMem)
      end if

c   o assign iUHF
      if (iflags(11).eq.0) then
         iUHF = 0
      else
         iUHF = 1
      end if

c   o initialize the chemical system (held in common blocks)
      call aces_init_chemsys


      return
c     end subroutine aces_init
      end

c ----------------------------------------------------------------------
      blockdata aces_bd_aces_reflag
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/
      data                 xFlag /0/
      end
c ----------------------------------------------------------------------

