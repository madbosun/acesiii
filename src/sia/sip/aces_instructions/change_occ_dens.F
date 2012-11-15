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
      subroutine change_occ_dens(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c
c     Changes the occupation in the density array
c
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
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
      integer*8 indblk, get_block_index
      integer stack
      
      integer comm

      integer val1(mx_array_index), val2(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer*8 addr, get_index_from_base

      integer sind,evec_array,evec_type
      double precision sval

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
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      array = op(c_result_array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         index = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg = index_table(c_current_seg,index)

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

c---------------------------------------------------------------------------
c   Get array data address.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,array) .eq. static_array) then
         addr = address_table(array)
         indblk = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack = array_table(c_array_stack,array)
         indblk = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif 

      if (array_table(c_array_type,array) .eq. static_array) then
         na1 = array_table(c_index_range1, array)
         na2 = array_table(c_index_range2,array)
         ni1 = array_table(c_index_range1+1,array)
         ni2 = array_table(c_index_range2+1,array)
      else
         na1 = val1(1)
         na2 = val2(1)
         ni1 = val1(2)
         ni2 = val2(2)
      endif 

c---------------------------------------------------------------------------
c   LOCATE THE SCALAR (ie. Orbital to change)
c---------------------------------------------------------------------------

      evec_array = op(c_op1_array)
      evec_type  = array_table(c_array_type, evec_array)
      if ( evec_type .ne. scalar_value) then
         print *,'Error: The sind argument must be a scalar.'
         call abort_job()
      endif

      sind = array_table(c_scalar_index, evec_array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
        print *,'Scalar table index out of range in change_occ,',
     *           'line ',current_line
        print *,'Index for array ',evec_array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif
      sval = scalar_table(sind)


      call change_the_dens_occ(x(indblk), val1(1),val2(1),type(1),
     *                         na1, na2,
     *                         val1(2),val2(2),type(2),ni1,ni2,
     *                         sval)
      
      return
      end
C------------------------------------------------------------------------C
C                                                                        C
C                                                                        C
C                                                                        C
C------------------------------------------------------------------------C

      subroutine change_the_dens_occ(x, p1,p2,ptype,np1,np2,
     *                       q1,q2,qtype,nq1,nq2,orbval)
      implicit none
      include 'interpreter.h'

      integer p1,p2,q1,q2
      integer ptype, qtype
      integer np1, np2, nq1, nq2
      double precision orbval
      double precision x(np1:np2,nq1:nq2)

      integer p, q, iorbval

      iorbval = INT(orbval)
      DO p = p1, p2
      DO q = q1, q2
         IF ( (p .EQ. iorbval) .AND. (p .EQ. q) ) x (p,q) = 0.0D0
      END DO
      END DO

      return
      end
