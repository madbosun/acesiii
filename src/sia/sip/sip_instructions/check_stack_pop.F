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
     *                           nindex_table, block_map_table) 
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
      include 'timerz.h' 

      integer array_table, narray_table, index_table, nindex_table 
      integer block_map_table(*) 

      integer array, block, type, request, instruction_timer,
     *        comm_timer, iblkndx, n_used, n_free, min_blocks  
      integer stack, ncount, min_block, max_siter 
      integer ierr
      integer status(MPI_STATUS_SIZE)
      logical flag

c     if (me .eq. 0) then 
c     write(6,*) 'Checking stack population ' 
c     write(6,*) 'Number of stacks:', nblkmgr_stacks  
c     endif 

      min_blocks = 7 
      max_siter  = 5000001  

      nblocks = 0 
      do i = nblkmgr_stacks, 1, -1 ! nblkmgr_stacks-1, -1 
            stack = i 
            nblocks = nblocks + nblocks_stack(stack) 
            ncount = 0 
            min_block = min(min_blocks,nblocks_stack(stack))  
c           min_block = min_block/3    
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
10             continue 
               ncount = ncount + 1 

               call scrub_from_stack(stack, array_table,
     *              narray_table, index_table, nindex_table,
     *              block_map_table, ierr)

               if (ierr .eq. 0) n_free = n_free + 1  
               if (n_free .ge. min_block) go to 11  

               call reclaim_persistent_block_from_stack(stack,
     *              array_table,narray_table, index_table, nindex_table,
     *              block_map_table, ierr)

               if (ierr .eq. 0) n_free = n_free + 1  
               if (n_free .ge. min_block) go to 11  

               if (ncount .lt. max_siter) go to 10 
            endif  

c     write(6,*) ' ME STACK NFREE :', me, i, n_free, 'AFTER :', ncount  
            write(6,*) ' nused nfree ', i, n_used, n_free, 
     *                 ' ntot ', nblocks_stack(i), 
     *                 ' start ', stack_start(i),    
     *                 ' iblkndx', iblkndx , 
     *                 ' after ', ncount, 'iterations on proc', me   
11       continue 
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
      integer array, type 
      integer n_allocate

c----------------------------------------------------------------------- 
c     Initialize the stack population 
c----------------------------------------------------------------------- 

      do i = 1, nstacks
         need_stack(i) = 0
      enddo
      n_allocate = 0

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

         do j = 1, 171 
            stack  = array_table(c_array_stack, j)
            write(6,*) j, stack 
         enddo 
         STOP 

c----------------------------------------------------------------------- 
c        Check for allocation of local array  
c----------------------------------------------------------------------- 

         if (opcode .eq. allocate_op) then
            array  = optable(c_result_array, i)
            type   = array_table(c_array_type,array) 
            if (type .ne. local_array) then
                print *,'Error: Deallocate instruction requires 
     *                   a local_array'
                print *,'Array, type = ',array,type
            endif 
            stack  = array_table(c_array_stack, array)
            nind   = array_table(c_nindex, array)
            if(array.eq.171) write(6,*) 'ARRAY NIND STACK:', array, 
     *        nind, stack 
            nblock = 1
            nwild  = 0
            do j = 1, nind
               ind = array_table(c_index_array1+j-1, array)
               nseg = index_table(c_nsegments, ind)
               if (optable(c_ind1+j-1,i) .eq. wildcard_indicator) then
                   nblock = nblock*nseg
                   nwild = nwild + 1 
               endif
            enddo
            need_stack(stack) = need_stack(stack) + nblock
            n_allocate = n_allocate + 1
         endif ! allocate 

c----------------------------------------------------------------------- 
c        Check for array assignment  
c----------------------------------------------------------------------- 

         if (opcode .eq. assignment_op) then
            array               = optable(c_result_array, i)
            type                = array_table(c_array_type,array) 
            if (type .ne. scalar_value) then 
                stack               = array_table(c_array_stack, array)
                need_stack(stack)   = need_stack(stack) + 3
                need_stack(nstacks) = need_stack(nstacks) + 3
            endif 
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

c----------------------------------------------------------------------- 
c        Check for deallocation of local array  
c----------------------------------------------------------------------- 

         if (opcode .eq. deallocate_op) then
            n_allocate = n_allocate - 1
         endif ! deallocate 

      enddo ! i 

100   continue

      current_line = optable(c_lineno,iopsave) 
      if (n_allocate .ne. 0) then
         write(6,*) ' You have not deallocated all local arrays in 
     *                the pardo at line :', current_line
      endif

      write(6,*) ' Stack requirements for pardo at line', current_line
      do i = 1, nstacks
         write(6,*) '    ', i, need_stack(i) 
      enddo

      return
      end

