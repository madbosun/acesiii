c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It defines if current point g belongs to certain processor   |     
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine grid_point_range(array_table,
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
      include 'pst_functions.h'
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

      integer i, j, k,qq,qq1,is_it_worker(max_procsx)
      integer nind, nsend,ind1,w_worker(max_procsx)
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
        
c-------------------------------------------------------------------

c--------------------------------------------------------------------

c------------Information about scalar--------------------------------------

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      g=scalar_table(sind) 

      our_workers = 0 
      do i = 1, nprocs
         if (pst_get_company(i-1) .eq. io_company_id) then
c           write(6,*) ' IO SERVER :', i  
            is_it_worker(i) = 0
         else
c           write(6,*) ' WORKER :', i  
            is_it_worker(i) = 1
            our_workers = our_workers + 1
            w_worker(i) = our_workers
         endif
      enddo


      ind1 = 1
      qq = (g-1)/our_workers
      qq1 = g - qq*our_workers
c---------------------------------------------------------------------------
c   Compute the block of integrals.
c---------------------------------------------------------------------------
c      write(6,*) "Current grid poit",g,'<--',me,qq1
c------------Information about scalar--------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      sind =  array_table(c_scalar_index, array)

      if(qq1.eq.w_worker(me+1)) ind1 = 1
      if(qq1.ne.w_worker(me+1)) ind1 = 0 

      dind1 = ind1
      scalar_table(sind) = dind1

c      write(6,*) '-->',g, w_worker(me+1), ind1,"<-",me+1 
      

      return
      end

