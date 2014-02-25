c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It copies eigenvectors of small Davidson matrix vectors      |
c  from block to static arrays DAV_R_AUX and DAV_L_AUX          | 
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine setauxvec (array_table,
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
      include 'daveom.h'
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
      integer*8 integral_scr,daddr,daddr1,iarray1
      integer seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag,ind(mx_array_index)
      integer stack,g,n1,n2,m1,m2

      double precision x(1),int_type_d,dind1
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

c-------------------DAV_R_aux array---------------------------------

      array      = op(c_result_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)


      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index
c-------------------DAV_L_aux array---------------------------------
      array = op(c_op1_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr1 = address_table(array)
      iarray1  = get_index_from_base(daddr1, x, 2)

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      m1 = index_table(c_index_size, ind(1))  ! pick up length of index
      m2 = index_table(c_index_size, ind(2))  ! pick up length of index
c---------------------------------------------------------------------
      call setauxvectors(x(iarray),n1,n2,x(iarray1),m1,m2 ) 
c---------------------------------------------------------------------


      return
      end


c#####################################################################
c#####################################################################
      subroutine setauxvectors(DAV_R_aux,n1,n2,DAV_L_aux,m1,m2) 
      IMPLICIT NONE
      include 'parallel_info.h'
      include 'daveom.h'

      integer i,j,n1,n2,m1,m2,izz,k
      double precision  ppp,
     & DAV_R_aux(n2,n2),DAV_L_aux(m2,m2) 

      izz = nroots*iter + nroots   

      do i=1,izz 
       do j=1,izz
       DAV_R_aux(i,j) = DavRvec(i,j)            
       DAV_L_aux(i,j) = DavLvec(i,j)
       enddo
      enddo

c--------Test-----------------------------

c      do i=1,izz  
c      do j=1,izz 
c         ppp = 0.d0
c         do k=1,izz
c         ppp = ppp+DavLvec(k,i)*DavRvec(k,j)  
c         enddo 
c         if(me.eq.0) write(6,*) ppp,'<--',i,j  
c      enddo    
c      enddo
c-----------------------------------------
      return
      end   

