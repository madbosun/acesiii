
c This routine defines the default values for the common blocks in
c cache.com.

      blockdata aces_bd_cache
      implicit none
c cache.com : begin

c These common blocks contain global information about the automatic file cache.
c getlst and putlst REQUIRE a cache, hence the term 'automatic' (compared to the
c auxiliary cache controlled by /auxcache/quikget).

c#define _CACHE_BYPASS /* bypasses cache on reading/writing of full records */
c#define _CACHE_HIST
c#define _CACHE_HIST_VERBOSE


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
      data cachemiss, cacheskip,
     &     cacheread, cachewrite, cachewriteback /5*0/
      data bCacheUp /.false./
      end

