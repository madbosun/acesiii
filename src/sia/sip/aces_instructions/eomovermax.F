c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs selection of maximal overlap from the overlap    |
C      table of current eom root                                |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine eomovermax
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
      include 'epsilon.h'
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
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, k
      integer array, index, nindex, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer*8 indblk1, indblk2, get_block_index,iarray,daddr
      integer stack
      integer nsend
      integer array_type, evec_type
      integer comm
      integer sind, ind(mx_array_index),root_index
      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2,n1,n2
      integer*8 addr, get_index_from_base
      double precision x(1)


      double precision overmax_index

#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------|
c   Locate the data for the input matrix which is to be analyses.           |
c   First array: set of overlaps for current root                           |
c---------------------------------------------------------------------------|

      array      = op(c_result_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      root_index = scalar_table(sind)  


c---------------------------------------------------------------------------|
c               Selection of maximal root                                   |
c---------------------------------------------------------------------------|

                 IF(RL.eq.1) then
                 overmax_index = DavOverRmax(root_index)
                 ENDIF

                 IF(RL.eq.0) then
                 overmax_index = DavOverLmax(root_index)
                 ENDIF

        

c---------------------------------------------------------------------------|
c   Scalar part: the output scalar will contain index of root               |
c                   with maximal overlap                                    |
c---------------------------------------------------------------------------|

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      scalar_table(sind) = overmax_index


      return
      end
c===========================================================================|
c===========================================================================|



C---------------------------------------------------------------------------|





