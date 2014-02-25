c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs calculation of values of molecular orbitals      |
C  z-derivatives at current grid point g                        |     
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine f12_2center_mzobj(array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)


c-------------------------------------------------------------------------- 
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'F12.h'
      include 'F12_gr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

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

      integer i, j, k
      integer nind, nsend,n_orb
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 iarray, scfevec, indblk, get_block_index
      integer*8 integral_scr
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack

      double precision x(1),int_type_p
      double precision dbuf(1), temp,g
      integer ibuf(1)
      integer int_type, obj_type,g1
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
      direct_flag = .false.

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c-------------------------------------------------------------------
        
c-------------------------------------------------------------------

c--------------------------------------------------------------------
c  At this step we create array aoint if it was not created before  |
c--------------------------------------------------------------------

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

      stack = array_table(c_array_stack,array)
      indblk = get_block_index(array, blk, stack,
     *                         blkndx, x, .true.)

c--------------------------------------------------------------------------
c   Get the segment ranges for each index of the integral array.
c--------------------------------------------------------------------------

      block = find_current_block(array, array_table(1,array),
     *                          index_table, nindex_table,
     *                          segment_table, nsegment_table,
     *                          block_map_table, blkndx)

      stack = array_table(c_array_stack, array)
      indblk = get_block_data_index(array, block, stack,
     *                              blkndx, x)




c---------------------------------------------------------------------------
c--------------------------------------------------------------------------
c   Get the segment ranges for each index of the integral array.
c--------------------------------------------------------------------------

      nind = array_table(c_nindex,array)
      if (nind .ne. 1) then
         print *,'Error: must be only 1 index. This sucks ',
     *        array,' which has ',nind,' indices.'
         print *,'       This routine requires 1 indices.'
         call abort_job()
      endif

      do i = 1, nind
         ind = array_table(c_index_array1+i-1,array)
         seg = index_table(c_current_seg,ind)
         call get_index_segment(ind, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

c------------Information about scalar--------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. scalar_value) then
         print *,'Error: The Second argument in get_ijk  
     *            must be a scalar.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in get_ijk, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      g = scalar_table(sind)
      g1 = g
c---------------------------------------------------------------------------

c           write(6,*) val1(1),val2(1),'<----'

      array = 2 ! SCF vectors alpha   
      array_type = array_table(c_array_type, array)

      if (array_type .ne. static_array) then
         print *,'Error: SCF vectors arrays in integral calc
     *            must be static.'
         call abort_job()
      endif

      nindex      = array_table(c_nindex, array)

      if (nindex .ne. 2) then
         print *,'Error: SCF vectors array must have ',
     *           'two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         indb(i) = index_table(c_index_size, na)  ! pick up length
      enddo

      scfevec = get_index_from_base(address_table(array), x, 2)
      n_orb = indb(1)
C----------------------------------------------------------------------|
c       Calculation of values of molecular orbitals z-derivatives at   |               
c   of given segment  at current point g                               |
                 obj_type = 4
c-----------------------------------------------------------------------
              
      if (.not. dryrun) then 
       call twocenter_obj(n_orb,x(scfevec),g1,val1(1),val2(1),x(indblk),
     & obj_type)
         
       endif


      return
      end

