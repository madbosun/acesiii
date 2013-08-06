
c This routine retrieves physical records from list files.
c It is a primitive of getlst and should not be called directly.

      subroutine getlst_io(ZList,iFile,iRecNum,iWrdNdx,nWords)
      implicit none

c ARGUMENTS
      integer ZList(*), iFile, iRecNum, iWrdNdx, nWords

c EXTERNAL FUNCTIONS
      integer iieq, iiamin, iiamax

c INTERNAL VARIABLES
      integer iFileOut, iRecOut
      integer iTmp, iLoc, iNdx, iVictim
      logical bNotDone

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

c ----------------------------------------------------------------------

      iTmp = 0
c   o assert cache subsystem is up
      if (.not.bCacheUp) then
         print *, '@GETLST_IO: Assertion failed.'
         print *, '   bCacheUp = ',bCacheUp
         iTmp = 1
      end if
c   o assert iFile, iRecNum, and iWrdNdx are all natural
      if ((iFile.lt.1).or.(iRecNum.lt.1).or.(iWrdNdx.lt.1)) then
         print *, '@GETLST_IO: Assertion failed.'
         print *, '   iFile   = ',iFile
         print *, '   iRecNum = ',iRecNum
         print *, '   iWrdNdx = ',iWrdNdx
         iTmp = 1
      end if
c   o assert nWords >= 0
      if (nWords.lt.0) then
         print *, '@GETLST_IO: Assertion failed.'
         print *, '   nWords = ',nWords
         iTmp = 1
      end if
c   o assert maximum element is in physical record
      if ((iWrdNdx-1+nWords).gt.iprcwd) then
         print *, '@GETLST_IO: Assertion failed.'
         print *, '   iWrdNdx = ',iWrdNdx
         print *, '   nWords  = ',nWords
         print *, '   iprcwd  = ',iprcwd
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)

      if ((iWrdNdx.lt.1).or.(nWords.lt.1)) return

c ----------------------------------------------------------------------

c   o increment the cachetime
      cachetime = cachetime + 1
      if (cachetime.lt.0) then
c      o scale-back the lrustats
         iLoc = iiamin(cachnum,lrustats,1)
         iTmp = lrustats(iLoc)
         do iLoc = 1, cachnum
            lrustats(iLoc) = lrustats(iLoc) - iTmp
         end do
c      o set cachetime to the next integer
         iLoc = iiamax(cachnum,lrustats,1)
         cachetime = lrustats(iLoc) + 1
      end if

c   o check to see if the file/record is in the cache
      bNotDone = .true.
      iLoc = cachnum
      do while (bNotDone.and.(iLoc.ne.0))
         if (cachrec(iLoc).ne.iRecNum) then
            iLoc = iLoc - 1
         else
            if (cachfil(iLoc).eq.iFile) then
               bNotDone = .false.
            else
               iLoc = iLoc - 1
            end if
         end if
      end do

c   o If the correct record is in the cache, then copy it to ZList.
      if (iLoc.ne.0) then
         iNdx = cachndx(iLoc)-1+iWrdNdx
         call icopy(nWords,icache(iNdx),1,ZList,1)
         lrustats(iLoc) = cachetime
      else


c   o the record must be picked up from disk
c     Dump the least recently used record to disk if it has been modified;
c     otherwise, just copy over it.
      iVictim = iiamin(cachnum,lrustats,1)
      iNdx = cachndx(iVictim)
      if (cachmod(iVictim).ne.0) then
         iFileOut = cachfil(iVictim)
         iRecOut  = cachrec(iVictim)
         call aces_io_write(iFileOut,iRecOut,icache(iNdx),iprcwd)
         cachmod(iVictim) = 0
      end if
      call aces_io_read(iFile,iRecNum,icache(iNdx),iprcwd)
      iNdx = iNdx-1+iWrdNdx
      call icopy(nWords,icache(iNdx),1,ZList,1)
      cachfil(iVictim) = iFile
      cachrec(iVictim) = iRecNum
      lrustats(iVictim) = cachetime


c     end if (iLoc.ne.0)
      end if

      return
c     end subroutine getlst_io
      end

