C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine index_match(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c--------------------------------------------------------------------------
c  In this routine the first array is copied into the second.  
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
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
      integer nscalar_table, array_type
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, k, sind 
      integer array, index, nindex, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer*8 indblk1, indblk2, get_block_index
      integer stack
      integer nsend 

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index) 
      integer type(mx_array_index) 
      integer na1, na2, ni1, ni2
      integer*8 addr, get_index_from_base
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Detrmine the parameters of the first array: c_result_array 
c---------------------------------------------------------------------------

      array = op(c_result_array)
      nindex = array_table(c_nindex, array) 

      if (nindex .ne. 2) then 
         print *,'ERROR: The array in index match MUST have two 
     *            indeces, yours has', nindex
         call abort_job() 
      endif 

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         sindex(i) = index 
         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, fop1(i), fop2(i))
      enddo


      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk1 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk1 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif

c---------------------------------------------------------------------------
c   Detrmine the parameters of the second array: c_op1_array 
c---------------------------------------------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. scalar_value) then
         print *,'Error: Second array in index_match must be a scalar.'
         call abort_job()
      endif

      sind = array_table(c_scalar_index, array)
  
      scalar_table(sind) = 0.0 
      do i = fop1(2), fop2(2) 
         if (i .eq. fop1(1)) scalar_table(sind) = 1.0 
      enddo 

c     write(6,*) 'INDEX_MATCH:', fop1(1), fop1(2), fop1(2), fop2(2), 
c    *            scalar_table(sind)  

      return
      end

