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
      subroutine retrieve_block_info(x, array, 
     *                      array_table, narray_table, 
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op,
     *                      indblk, nindex, type,
     *                      bval, eval, bdim, edim)
c--------------------------------------------------------------------------
c   Determines necessary data for the current block of a given array.
c
c   x              Base address used to index the block's address.
c   array          Array being referenced.
c   array_table    Table of array information
c   narray_table   Number of array table entries.
c   index_table    Index table.
c   nindex_table   Number of index table entries.
c   segment_table  Segment table.
c   nsegment_table Number of segment table entries.
c   block_map_table Block map table
c   nblock_map_table Number of block_map_table entries
c   scalar_table   Table of scalar values.
c   nscalar_table  Number of entries in scalar_table.
c   address_table  Table of addresses of static arrays. (integer*8)
c   op             Current instruction.
c   indblk         Return index of the block, relative to x. (integter*8)
c   nindex         Return number of indices in the array.
c   type           Return array of iindex types.
c   bval           Return array of beginning index values.
c   eval           Return array of ending index values.
c   bdim           Return array of beginning index dimensions.
c   edim           Return array of ending index dimensions.
c
c   If the array is a static array, the bdim/edim values will be different from 
c   the bval/eval values.  Otherwise, they will be the same.
c
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'

      integer narray_table, nindex_table, nsegment_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry,*)
      integer nblock_map_table
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(*)
      integer op(loptable_entry)

      integer array
      integer nindex
      integer*8 indblk
      integer bval(mx_array_index), eval(mx_array_index) 
      integer bdim(mx_array_index), edim(mx_array_index) 
      integer type(mx_array_index)

      integer i, stack, block, index, seg, blkndx, ierr, blk, flag
      integer find_current_block
      integer array_type, ind
      integer*8 addr
      integer*8 get_block_index, get_index_from_base
      integer*8 c_loc64
      integer*8 ixx
      logical direct_flag

      double precision x(1)

c---------------------------------------------------------------------------
c   Check for special case of a scalar variable.
c---------------------------------------------------------------------------

      array_type = array_table(c_array_type, array)
      if (array_type .eq. scalar_value) then
         ind =  array_table(c_scalar_index, array)
         if (ind .lt. 1 .or. ind .gt. nscalar_table) then
            print *,'Scalar table index out of range in print_scalar, ',
     *              'line ',current_line
            print *,'Index for array ',array,' is ',ind,' should be ',
     *              'between 1 and ',nscalar_table
            call abort_job()
         endif

c---------------------------------------------------------------------------
c   Calculate the index of the scalar value relative to base address "x".
c---------------------------------------------------------------------------

         ixx = 1
         addr = c_loc64(scalar_table(ind), ixx, 1)   
         indblk = get_index_from_base(addr, x, 2)  ! double precision index.

c----------------------------------------------------------------------------
c   Return nindex = 0 as a signal to the instructiun that this is a scalar.
c----------------------------------------------------------------------------

         nindex = 0
         return
      endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         index = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg = index_table(c_current_seg,index)
         if (seg .eq. undefined_segment) then
            if (array_table(c_array_type,array) .eq. 
     *                                   static_array) then

c-------------------------------------------------------------------------
c   Static array, using entire array as argument.  Set range of values to 
c   the entire range of the index, rather than one particular segment.
c-------------------------------------------------------------------------

               bval(i) = array_table(c_index_range1+i-1, array)
               eval(i) = array_table(c_index_range2+i-1, array)
            else

c--------------------------------------------------------------------------
c   Non-static array with undefined current segment.  Abort the job.
c--------------------------------------------------------------------------

              print *,'Array ',array,' is used with undefined segment',
     *                ' for index ',index
              call abort_job() 
            endif
         else

c-------------------------------------------------------------------------
c   Normal usage.  Get segment ranges, preset the dimensional ranges.
c-------------------------------------------------------------------------

            call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, bval(i), eval(i))
            bdim(i) = bval(i)
            edim(i) = eval(i)
         endif   
      enddo

c---------------------------------------------------------------------------
c   Get array data address.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,array) .eq. static_array) then
         addr = address_table(array)
         indblk = get_index_from_base(addr, x, 2)

c----------------------------------------------------------------------------
c   Override the dimension fl ranges for each index from array_table info.
c----------------------------------------------------------------------------

         do i = 1, nindex
            bdim(i) = array_table(c_index_range1+i-1, array)
            edim(i) = array_table(c_index_range2+i-1,array)
         enddo
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)
    
         if (blkndx .eq. 0) then
c----------------------------------------------------------------------------
c   Block does not already exist.  We must create it from scratch (i. e., 
c   bring the block into scope).
c----------------------------------------------------------------------------

            direct_flag = .true.
            call create_current_block(array,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., direct_flag, blk, ierr)
            blkndx = ierr
            call get_block_computed_flag(array, blk, blkndx, flag)
            if (flag .eq. 0) then
               call set_opblock(array, blk, blkndx, op)
               call set_block_computed_flag(array, blk, blkndx, 1)
            endif
         endif


         stack = array_table(c_array_stack,array)
         indblk = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif 

      return
      end

