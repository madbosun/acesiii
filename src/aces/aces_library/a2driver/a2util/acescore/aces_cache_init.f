
c This routine reserves space at iCore(i0) for an I/O cache.

c The address of iCore(i0) is saved before returning; therefore, there is no
c significance placed on these two variables besides pointing to the usable
c memory heap.

c INPUT
c int ICORE(*) : the memory anchor used to place the cache

c INTPUT/OUTPUT
c int I0 : on input, this is the ICORE index that places the cache
c        : on output, this points to the first integer after the cache
c int IMEM : the amount of usable memory at ICORE(I0)
c int NPRCWD : the desired/actual integer-length of a physical record
c int NSLOTS : the desired/actual number of cache slots








c#define _DEBUG_ACES_CACHE_INIT

      subroutine aces_cache_init(iCore,i0,iMem,nPRcWd,nSlots)
      implicit none

c ARGUMENTS
      integer iCore(*), i0, iMem, nPRcWd, nSlots

c EXTERNAL FUNCTIONS
      INTEGER c_adr
      external c_adr

c INTERNAL VARIABLES
      double precision dTmp
      INTEGER zTmp
      integer i, ndx, iTmp

c COMMON BLOCKS
c cache.com : begin

c These common blocks contain global information about the automatic file cache.
c getlst and putlst REQUIRE a cache, hence the term 'automatic' (compared to the
c auxiliary cache controlled by /auxcache/quikget).

c#define _CACHE_BYPASS /* bypasses cache on reading/writing of full records */
c#define _CACHE_HIST
c#define _CACHE_HIST_VERBOSE

      external aces_bd_cache

c icache     : the anchor used to address each cache slot
c cachnum    : the number of usable cache slots
c cachrec(i) : the index of the physical record cached by the data in slot i
c cachfil(i) : the external file unit number that stores the data in slot i
c cachndx(i) : the icache index of slot i
c cachmod(i) : a modification flag used to trigger a writeback

c cachetime   : a cache-event counter
c lrustats(i) : the last 'time' slot i was accessed

      integer        icache(1), cachnum,
     &               cachrec(128),
     &               cachfil(128),
     &               cachndx(128),
     &               cachmod(128)
      common /cache/ icache, cachnum,
     &               cachrec,
     &               cachfil,
     &               cachndx,
     &               cachmod
      save   /cache/

      integer           cachetime, lrustats(128)
      common /cachelru/ cachetime, lrustats
      save   /cachelru/

c cachemiss      : measures cache misses
c cacheskip      : measures read and write bypasses
c cacheread      : measures read hits
c cachewrite     : measures write hits
c cachewriteback : measures writes-back of dirty slots

      integer             cachemiss, cacheskip,
     &                    cacheread, cachewrite, cachewriteback
      common /cache_hist/ cachemiss, cacheskip,
     &                    cacheread, cachewrite, cachewriteback
      save   /cache_hist/

c bCacheUp : a flag for bombing in get/putlst if there is no I/O cache

      logical              bCacheUp
      common /cache_flags/ bCacheUp
      save   /cache_flags/

c cache.com : end
c filspc.com : begin

c This common block contains the dimensions of the physical records used by the
c MOIO storage files.

c iprcln : the    byte-length of a physical record
c iprcwd : the integer-length of a physical record
c iprcrl : the    recl-length of a physical record

      integer         iprcln, iprcwd, iprcrl
      common /filspc/ iprcln, iprcwd, iprcrl
      save   /filspc/

c filspc.com : end


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



c icdacc.com : begin
c Nevin 8-30-95 added record length common to facilitate change from
c Bytes to Words for SGI and DecAlpha
      integer         idaccm
      common /icdacc/ idaccm
c icdacc.com : end

c ----------------------------------------------------------------------

c   o the cache subsystem never turns off
      if (bCacheUp) return

c   o turn on the cache subsystem flag
      bCacheUp = .true.

c   o make sure iMem is natural
      if (iMem.lt.1) then
         print *, '@ACES_CACHE_INIT: There is no place for the cache.'
         print *, '   iMem = ',iMem
         call aces_exit(1)
      end if

c   o initialize the filspc common block
c     (set the length of a physical record and cache slot)
      if (nPRcWd.lt.1) then
         print *, '@ACES_CACHE_INIT: The length of a cache slot must ',
     &            'be at least 1 integer long.'
         print *, '                  requested length = ',nPRcWd
         call aces_exit(1)
      end if
      iprcln = nPRcWd*iintln
      iprcwd = nPRcWd
      iprcrl = nPRcWd*idaccm
      if (iand(iprcln,ifltln-1).ne.0) then
         print '(/)'
         print *, '@ACES_CACHE_INIT: WARNING - The length of a cache ',
     &            'slot is not a whole multiple'
         print *, '                  of floating-point numbers. This ',
     &            'will negatively impact'
         print *, '                  performance.'
         print '(/)'
      end if

c   o set the number of cache slots
      if (nSlots.lt.1) then
         print *, '@ACES_CACHE_INIT: There must be at least 1 cache ',
     &            'slot.'
         print *, '                  requested number = ',nSlots
         call aces_exit(1)
      end if
      if (nSlots.gt.128) then
         print *, '@ACES_CACHE_INIT: There can be at most ',
     &            128,' cache slots.'
         print *, '                  requested number = ',nSlots
         print *, '                  resetting to ',128
         nSlots = 128
      end if
      cachnum = nSlots

c   o use ndx to point to the first double boundary in iCore(i0)
      zTmp = c_adr(iCore(i0))
      if (iand(zTmp,ifltln-1).ne.0) then
         i0 = i0 + 1
         zTmp = c_adr(iCore(i0))
         if (iand(zTmp,ifltln-1).ne.0) then
            print *, '@ACES_CACHE_INIT: The core memory is not aligned.'
            print *, '   &icore[0] = ',c_adr(iCore(1))
            call aces_exit(1)
         end if
         iMem = iMem - 1
      end if
      zTmp = zTmp - c_adr(icache(1))
      if (iand(zTmp,iintln-1).ne.0) then
         print *, '@ACES_CACHE_INIT: The cache anchor is not aligned.'
         print *, '   &icache[0] = ',c_adr(icache(1))
         call aces_exit(1)
      end if
      zTmp = 1+(zTmp/iintln)
      ndx = zTmp
      if (ndx.ne.zTmp) then
         print *, '@ACES_CACHE_INIT: Assertion failed.'
         print *, '                  The I/O cache is not addressable.'
         call aces_exit(1)
      end if
      do i = 1, cachnum
         cachndx(i) = ndx
         ndx = ndx + iprcwd
      end do

c   o move i0 to a double boundary and decrease iMem
      i = iprcwd * cachnum
      i0 = i0 + i
      zTmp = c_adr(iCore(i0))
      if (iand(zTmp,ifltln-1).ne.0) then
         i0 = i0 + 1
         iMem = iMem - 1
      end if
      iMem = iMem - i
      if (iMem.lt.1) then
         print *, '@ACES_CACHE_INIT: The total cache size must be less',
     &            ' than the total core size.'
         print *, '                  cache slots = ',cachnum
         print *, '                  slot size   = ',iprcwd,' integers'
         print *, '                  core size   = ',iMem,' integers'
         call aces_exit(1)
      end if

c   o (re)set the cache statistics
      call aces_cache_reset


      return
c     end subroutine aces_cache_init
      end

