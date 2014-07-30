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
      subroutine frag_index_range(array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)


c-------------------------------------------------------------------------- 


c--------------------------------------------------------------------------
c Input: x1=which fragment
c        x2=type
c           1=aoindex 
c           2=moindex 
c           3=moindex occupied  
c           4=moindex virtual   
c           5=moaindex 
c           6=moaindex occupied  
c           7=moaindex virtual   
c           8=mobindex 
c           9=mobindex occupied  
c           10=mobindex virtual   
c Output: Beginning segment, Ending segment 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'
      include 'frag_def.h' 

      double precision x1, x2
      integer istart, iend   

c---------------Variables--------------------------------------------

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, array
      integer op(loptable_entry),sind
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table), get_index_from_base
      integer*8 get_block_data_index
      integer find_current_block, block

      integer array_type, nindex, inda(2), indb(2)

      logical direct_flag

      integer i, j, k,qq,qq1
      integer nind, nsend,ind1
      integer ierr, handle,our_workers
      integer blk, blkndx, maxblk
      integer*8 iarray, farray, indblk, get_block_index
      integer*8 integral_scr
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack,g

      double precision x(1),int_type_p,dind1
      double precision dbuf(1), temp
      integer ibuf(1)
      integer int_type
c----------------------------------------------------------------------
c----------------------------------------------------------------------
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (xptr, x)
      pointer (iptr, ibuf)
#else
      common dbuf, x
      equivalence (ibuf(1), dbuf(1))
#endif

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter

      array = op(c_result_array)

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c-------------------------------------------------------------------

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      x1 = scalar_table(sind)
C-------------------------------------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      x2 = scalar_table(sind)

C-------------------------------------------------------------------

      wfrag = x1 

      if (x2 .eq. 1.0) then 
         istart = 0 
         iend   = 0 
         if (wfrag .eq. 1) then 
            istart = 1 
            go to 10 
         endif 
         do i = 2, maxsegs  
            if ((w_frag_AO(i) .eq. wfrag) .and.
     *          (w_frag_AO(i-1) .ne. wfrag)) then 
                istart = i
                go to 10 
            endif 
         enddo 
10       continue 
         do i = istart, maxsegs  
            if ((w_frag_AO(i-1) .eq. wfrag) .and.
     *          (w_frag_AO(i) .ne. wfrag)) then 
                iend = i-1 
                go to 20 
            endif 
         enddo 
20       continue 
      endif 

      if (x2 .eq. 6.0) then 
         istart = 0 
         iend   = 0 
         if (wfrag .eq. 1) then 
            istart = 1 
            go to 11 
         endif 
         do i = 2, maxsegs  
            if ((w_frag_aocc(i) .eq. wfrag) .and.
     *          (w_frag_aocc(i-1) .ne. wfrag)) then 
                istart = i
                go to 11 
            endif 
         enddo 
11       continue 
         do i = istart, maxsegs  
            if ((w_frag_aocc(i-1) .eq. wfrag) .and.
     *          (w_frag_aocc(i) .ne. wfrag)) then 
                iend = i-1  
                go to 21 
            endif 
         enddo 
21       continue 
      endif 

      if (x2 .eq. 9.0) then 
         istart = 0 
         iend   = 0 
         if (wfrag .eq. 1) then 
            istart = 1 
            go to 12 
         endif 
         do i = 2, maxsegs  
            if ((w_frag_bocc(i) .eq. wfrag) .and.
     *          (w_frag_bocc(i-1) .ne. wfrag)) then 
                istart = i
                go to 12 
            endif 
         enddo 
12       continue 
         do i = istart, maxsegs  
            if ((w_frag_bocc(i-1) .eq. wfrag) .and.
     *          (w_frag_bocc(i) .ne. wfrag)) then 
                iend = i-1  
                go to 22 
            endif 
         enddo 
22       continue 
      endif 

      if (x2 .eq. 7.0) then 
         istart = 0 
         iend   = 0 
         if (wfrag .eq. 1) then 
            istart = 1 
            go to 13 
         endif 
         do i = 2, maxsegs  
            if ((w_frag_avirt(i) .eq. wfrag) .and.
     *          (w_frag_avirt(i-1) .ne. wfrag)) then 
                istart = i
                go to 13 
            endif 
         enddo 
13       continue 
         do i = istart, maxsegs  
            if ((w_frag_avirt(i-1) .eq. wfrag) .and.
     *          (w_frag_avirt(i) .ne. wfrag)) then 
                iend = i-1  
                go to 23 
            endif 
         enddo 
23       continue 
      endif 

      if (x2 .eq. 10.0) then 
         istart = 0 
         iend   = 0 
         if (wfrag .eq. 1) then 
            istart = 1 
            go to 14 
         endif 
         do i = 2, maxsegs  
            if ((w_frag_bvirt(i) .eq. wfrag) .and.
     *          (w_frag_bvirt(i-1) .ne. wfrag)) then 
                istart = i
                go to 14 
            endif 
         enddo 
14       continue 
         do i = istart, maxsegs  
            if ((w_frag_bvirt(i-1) .eq. wfrag) .and.
     *          (w_frag_bvirt(i) .ne. wfrag)) then 
                iend = i-1  
                go to 24 
            endif 
         enddo 
24       continue 
      endif 



      x1 = istart 
      x2 = iend 
c-------------------------------------------------------

c------------Information about scalar-------------------

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      scalar_table(sind) = x1
C-------------------------------------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      scalar_table(sind) = x2
c-------------------------------------------------------


c     write(6,*) ' AOSTART AOEND IN FRAG:', x1, x2 

      return
      end

