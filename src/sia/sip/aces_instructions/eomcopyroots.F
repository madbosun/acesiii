c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It copies roots from the coomon block to static array        |
c                in the EOM SIAL code                           |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine eomcopyroots
     +                 (array_table, 
     *                  narray_table, index_table,
     *                  nindex_table, segment_table, nsegment_table,
     *                  block_map_table, nblock_map_table,
     *                  scalar_table, nscalar_table, 
     *                  address_table, op)

c--------------------------------------------------------------------------
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'daveom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table 
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer type(mx_array_index) 
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer*8 iarray, ievec, get_index_from_base, daddr
      integer i, j, n1, n2, n3, n4  
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer sind, ind(mx_array_index)
      integer seg, index(mx_array_index), val1(mx_array_index), 
     *        val2(mx_array_index), junk

      integer block, blkndx
      integer find_current_block
      integer*8 get_block_index
      integer stack

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
#ifdef ALTIX
      dptr = dshptr
#endif




c----------------------------------------------------------------------------
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_result_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      iarray = get_index_from_base(address_table(array), x, 2)

      if (array .eq. 0) then
         print *,'Error: The first argument in return_sval  
     *            must be an array.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index


C----------------------------------------------------------------------------

      call eom_copyroot(x(iarray),1,n1,iter,nroots,DavRoots)

C----------------------------------------------------------------------------

      return 
      end 
C#############################################################################

      subroutine eom_copyroot(Dav_Roots,a1,a2,iter,nroots,DavRoots)
      IMPLICIT NONE

      integer i,j,a1,a2,b1,b2,iter,nroots,IWRK,INFO,LCI,KK,ICK
      double precision DAV_Roots(a1:a2)
      double precision DavRoots(20)
      include 'parallel_info.h'      

      do i=1,nroots
      DAV_Roots(i) = DavRoots(i)
      enddo

      return
      end
c-----------------------------------------------------------------------------
