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
      subroutine check_stack_pop(array_table, narray_table, index_table,
     *                           nindex_table, block_map_table, 
     *                           need_stack, nstacks) 
c-------------------------------------------------------------------------
c The population of the stacks is checked. If pop is too great a wait is
c enforced so the subsequent operations can proceed.  
c-------------------------------------------------------------------------

      implicit none
      integer i, j, nblocks, blkndx   

      include 'interpreter.h' 
      include 'blkmgr.h' 
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'server_monitor.h'
      include 'trace.h'

      integer array_table, narray_table, index_table, nindex_table 
      integer block_map_table(*) 

      integer array, block, type, request, instruction_timer,
     *        comm_timer, iblkndx, n_used, n_free, min_blocks  
      integer stack, ncount, min_block, max_siter 

      integer nstacks, need_stack(nstacks) 

      integer ierr
      integer status(MPI_STATUS_SIZE)
      logical flag

c     if (me .eq. 0) then 
c     write(6,*) 'Checking stack population ' 
c     write(6,*) 'Number of stacks:', nblkmgr_stacks  
c     endif 

      max_siter = 5000001  

      nblocks = 0 
      do i = nblkmgr_stacks, nblkmgr_stacks-1, -1 
            stack     = i 
            min_block = need_stack(stack)  
            nblocks   = nblocks + nblocks_stack(stack) 
            ncount = 0 
10          continue 
            ncount = ncount + 1 
            call find_free_stack(stack,iblkndx) 
            if (iblkndx .lt. 0) then 
               n_free = 0 
               n_used = nblocks_stack(stack) 
            else 
               n_free = iblkndx - stack_start(stack) + 1 
               n_used = nblocks_stack(stack) - n_free  
            endif 

            if (n_free .ge. min_block) go to 11  

            if (n_free .lt. min_block) then 
            call exec_thread_server(0) 

            call scrub_from_stack(stack, array_table,
     *           narray_table, index_table, nindex_table,
     *           block_map_table, ierr)

            call reclaim_persistent_block_from_stack(stack,
     *           array_table,narray_table, index_table, nindex_table,
     *           block_map_table, ierr)

            endif  

         if (ncount .lt. max_siter) go to 10 
         write(6,*) ' There are only',n_free,' free blocks on stack', 
     *                stack, 'on processor ', me 
11       continue 
c           write(6,*) ' nused nfree ', i, n_used, n_free, 
c    *                 ' ntot ', nblocks_stack(i), 
c    *                 ' start ', stack_start(i),    
c    *                 ' iblkndx', iblkndx , 
c    *                 ' after ', ncount, 'iterations'   
      enddo 

      return 
      end 

      subroutine find_min_stack_op_pardo(optable, noptable,
     *                   array_table,
     *                   narray_table,
     *                   array_labels, index_table, nindex_table,
     *                   segment_table, nsegment_table, block_map_table,
     *                   nblock_map_table,
     *                   scalar_table, nscalar_table, proctab,
     *                   address_table,
     *                   iopsave, need_stack, nstacks)
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_barrier_data.h'
      include 'scratchpad.h'
      include 'dbugcom.h'
      include 'where_table.h'
      include 'context.h'
      include 'checkpoint_data.h'
      include 'pst_functions.h'
      include 'int_gen_parms.h'
      include 'blkmgr.h'

      integer noptable, narray_table, nindex_table, nsegment_table
      integer nblock_map_table, nscalar_table
      integer comm
      integer optable(loptable_entry,noptable)
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer proctab(2,*)
      character*10 array_labels(narray_table)
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer ierr
      integer iopsave, iblk
      integer index, result_array, result_type
      integer block_map_entry(lblock_map_entry)
      integer op1_block_map_entry(lblock_map_entry)
      integer op2_block_map_entry(lblock_map_entry)
      integer opcode
      integer i, j, k, ind, nind, nseg  
      integer stack, nblock, nwild, nstacks, need_stack(nstacks)
      integer array 
      integer op(loptable_entry) 

c----------------------------------------------------------------------- 
c     Initialize the stack population 
c----------------------------------------------------------------------- 

      do i = 1, nstacks 
         need_stack(i) = 0 
      enddo 

c----------------------------------------------------------------------- 
c     Loop through optable between(opsave, end_pardo) counting blocks
c     needed.  
c----------------------------------------------------------------------- 

      opcode = optable(c_opcode, iopsave) 
      if (opcode .ne. pardo_op) then 
         write(6,*) ' Attemping to find stack population needed for a
     *                pardo but the opcode is:', opcode  
         call abort_job() 
      endif 

      do i = iopsave+1, noptable 

c----------------------------------------------------------------------- 
c        If end of pardo exit 
c----------------------------------------------------------------------- 

         opcode = optable(c_opcode, i) 
         if (opcode .eq. endpardo_op) go to 100 

c----------------------------------------------------------------------- 
c        Check for allocation of local array  
c----------------------------------------------------------------------- 

         if (opcode .eq. allocate_op) then 
            array  = optable(c_result_array, i) 
            nind   = array_table(c_array_type, array) 
            stack  = array_table(c_array_stack, array) 
            nblock = 1 
            nwild  = 0  
            do j = 1, nind 
               ind = array_table(c_index_array1+j-1, array) 
               nseg = index_table(c_nsegments, ind) 
               if (optable(c_ind1+j-1,i) .eq. wildcard_indicator) then
                   nblock = nblock*nseg
               endif 
            enddo 
            need_stack(stack) = need_stack(stack) + nblock 
         endif ! allocate 

c----------------------------------------------------------------------- 
c        Check for array assignment  
c----------------------------------------------------------------------- 

         if (opcode .eq. assignment_op) then 
            array               = optable(c_result_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 3 

            need_stack(nstacks) = need_stack(nstacks) + 3 
         endif ! assignment_op 

c----------------------------------------------------------------------- 
c        Check for contraction  
c----------------------------------------------------------------------- 

         if (opcode .eq. contraction_op) then 
            array               = optable(c_result_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 

            array               = optable(c_op1_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 

            array               = optable(c_op2_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 

            need_stack(nstacks) = need_stack(nstacks) + 3 
         endif ! assignment_op 

c----------------------------------------------------------------------- 
c        Check for summation  
c----------------------------------------------------------------------- 

         if ((opcode .eq. sum_op) .or. (opcode .eq. subtract_op)) then 
            array               = optable(c_result_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 

            array               = optable(c_op1_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 

            array               = optable(c_op2_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 

            need_stack(nstacks) = need_stack(nstacks) + 3 
         endif ! summation  

c----------------------------------------------------------------------- 
c        Check for get/request   
c----------------------------------------------------------------------- 

         if ((opcode .eq. get_op) .or. (opcode .eq. request_op)) then 
            array               = optable(c_result_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 
         endif ! get/request   

c----------------------------------------------------------------------- 
c        Check for put/prepare  
c----------------------------------------------------------------------- 

         if ((opcode .eq. put_op) .or. (opcode .eq. prepare_op)) then 
            array               = optable(c_result_array, i) 
            stack               = array_table(c_array_stack, array) 
            need_stack(stack)   = need_stack(stack) + 1 
         endif ! get/request   
 
      enddo ! i 

100   continue 

      return 
      end 
