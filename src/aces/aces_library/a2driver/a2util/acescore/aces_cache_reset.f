
c This routine resets the automatic cache data. It DOES NOT write-back
c any modified slots.

      subroutine aces_cache_reset
      implicit none

c INTERNAL VARIABLES
      integer i

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

c ----------------------------------------------------------------------

c   o make sure the cache subsystem is up
      if (.not.bCacheUp) return

c   o zero the cached addresses
      do i = 1, cachnum
         cachrec(i) = 0
      end do

c   o zero the LUNs
      do i = 1, cachnum
         cachfil(i) = 0
      end do

c   o reset the modification flags
      do i = 1, cachnum
         cachmod(i) = 0
      end do

c   o reset the cache event counter
      cachetime = 0

c   o zero the LRU statistics
      do i = 1, cachnum
         lrustats(i) = 0
      end do

      return
c     end subroutine aces_cache_reset
      end

