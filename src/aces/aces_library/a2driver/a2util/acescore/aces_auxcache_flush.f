
c This routine flushes any memory-resident list pointed to by quikget.
c The operation is also called a "write-back". This routine DOES NOT zero
c (destroy) the contents of quikget since doing so would perpetuate sloppy
c programming by losing the memory pointers.

      subroutine aces_auxcache_flush
      implicit none

c EXTERNAL FUNCTIONS
      integer aces_list_cols

c INTERNAL VARIABLES
      integer iFam, iGrp, iNdx, nCols

c COMMON BLOCKS


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end








c auxcache.com : begin

c The auxiliary cache is a programmer-controlled list cache. The dimensions are
c the same as those of the MOIO arrays. The programmer may load lists into icore
c memory and set quikget(?,?) to their icore addresses. When getlst and putlst
c operate on ANY list, quikget is checked to see if the list lives in icore.
c If so, the operation is performed on the in-core data instead of hitting the
c storage file(s).

c WARNING
c    There is no automatic updating of storage files with data in the auxiliary
c cache. If the memory-resident data is altered and must be stored to disk, then
c the quikget value must be destroyed (zeroed) and the actual location of the
c data must then be passed to putlst. The aces_auxcache_flush routine
c systematically stores the quikget values, calls putlst, and restores the
c quikget values.





      integer           quikget(10,500)
      common /auxcache/ quikget
      save   /auxcache/

c auxcache.com : end

c ----------------------------------------------------------------------

c   o write out any in-core list
      do iFam = 1, 500
      do iGrp = 1, 10
         if (quikget(iGrp,iFam).ne.0) then
            iNdx  = quikget(iGrp,iFam)
            nCols = aces_list_cols(iGrp,iFam)
            quikget(iGrp,iFam) = 0
            call putlst(icore(iNdx),1,nCols,1,iGrp,iFam)
            quikget(iGrp,iFam) = iNdx
         end if
      end do
      end do

      return
c     end subroutine aces_auxcache_flush
      end

