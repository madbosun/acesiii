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
      subroutine maxamp_print(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
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
      integer*8 indblk, get_block_index
      integer block, blkndx, seg
      integer find_current_block

      integer i, j, ii, source, company, comm, ierr
      integer source_type
      integer*8 isource 
      integer nindex_source, stack  
      integer n(mx_array_index),ind(mx_array_index) 
      integer sval1(mx_array_index),sval2(mx_array_index)  
      integer tval1(mx_array_index),tval2(mx_array_index)  
      integer msg(len_sip_server_message)
      integer status(MPI_STATUS_SIZE) 

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif


c--------------------------------------------------------------------------
c   Locate the data for both arrays.  
c--------------------------------------------------------------------------

      source = op(c_result_array)

c---------------------------------------------------------------------------
c   Look up source's address.
c---------------------------------------------------------------------------

      block  = find_current_block(source, array_table(1,source),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,source)
      isource = get_block_index(source, block, stack, blkndx, x, .true.) 

c--------------------------------------------------------------------------
c   Check the dimensions of both arrays both arrays.  
c--------------------------------------------------------------------------
      
      source_type = array_table(c_array_type, source)

      nindex_source = array_table(c_nindex, source)

      if (nindex_source .gt. 4) then 
         print *,'Arrays in maxamp_print must have no more than 
     *            4 indices'  
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source 
         ind(i) = array_table(c_index_array1+i-1,source)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, sval1(i), sval2(i))

      enddo


      if (nindex_source .eq. 4) then 
      call do_maprint4(x(isource),
     *                      sval1(1),sval2(1),sval1(2),sval2(2),
     *                      sval1(3),sval2(3),sval1(4),sval2(4)) 
      return 
      endif 

      if (nindex_source .eq. 2) then 
      call do_maprint2(x(isource),
     *                      sval1(1),sval2(1),sval1(2),sval2(2))
      return 
      endif 

      return
      end

      subroutine do_maprint4(x,i1,i2,j1,j2,k1,k2,l1,l2)
      implicit none 
      include 'interpreter.h'
      include 'parallel_info.h'

      integer i1,i2,j1,j2,k1,k2,l1,l2
      double precision x(i1:i2,j1:j2,k1:k2,l1:l2)

      integer a,b,c,d 
      double precision stemp

      open(66,file='summary.out')
c
       if (me .eq. 0)then
c
      do d = l1, l2
      do c = k1, k2
      do b = j1, j2
      do a = i1, i2
c
c     if(a .lt. c .and. b .lt. d)then
c
      stemp=dabs(x(a,b,c,d))
c
      if(stemp .gt. 0.1d-1)then
       write(66,1000) a, b, c, d, x(a,b,c,d)   
      endif
c
c     endif
c
      enddo
      enddo
      enddo
      enddo
c
       endif

 1000 FORMAT (4I8,D20.6)
c
      return 
      end 


      subroutine do_maprint2(x,i1,i2,j1,j2)
      implicit none 
      include 'interpreter.h'
      include 'parallel_info.h'

      integer i1,i2,j1,j2
      double precision x(i1:i2,j1:j2)

      integer a,b,c,d 
      double precision stemp

      open(66,file='summary.out')
c
       if (me .eq. 0)then
c
      do b = j1, j2
      do a = i1, i2
c
c     if(a .lt. c .and. b .lt. d)then
c
      stemp=dabs(x(a,b))
c
      if(stemp .gt. 0.1d-1)then
       write(66,1000) a, b, x(a,b)   
      endif
c
c     endif
c
      enddo
      enddo
c
       endif

 1000 FORMAT (2I8,D20.6)
c
      return 
      end 

