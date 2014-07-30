c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It makes diagonal elements of 3-body terms of EOM-IP matrix  |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|


      subroutine make_3body_diag_ip(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'trace.h'
      include 'epsilon.h'  
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

      integer i, j, ii, source, target, company, comm, ierr
      integer source_type, target_type
      integer nindex_source, nindex_target, stack
      integer*8 isource, itarget
      integer n(mx_array_index),ind(mx_array_index) 
      integer sval1(mx_array_index),sval2(mx_array_index)  
      integer tval1(mx_array_index),tval2(mx_array_index)  
      integer msg(len_sip_server_message)
      integer status(MPI_STATUS_SIZE) 
      integer type(mx_array_index)
   
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

      target = op(c_op1_array)
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

c---------------------------------------------------------------------------
c   Look up target's address.
c---------------------------------------------------------------------------

      block  = find_current_block(target, array_table(1,target),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,target)
      itarget = get_block_index(target, block, stack, blkndx, x, .true.) 

c--------------------------------------------------------------------------
c   Check the dimensions of both arrays both arrays.  
c--------------------------------------------------------------------------
      
      source_type = array_table(c_array_type, source)
      target_type = array_table(c_array_type, target)

      nindex_source = array_table(c_nindex, source)
      nindex_target = array_table(c_nindex, target)

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source 
         ind(i) = array_table(c_index_array1+i-1,source)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 
         type(i) = index_table(c_index_type,ind(i))

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, sval1(i), sval2(i))

      enddo

c-------------------------------------------------------------------------
c   Get segment ranges of the array target.
c-------------------------------------------------------------------------

      do i = 1, nindex_target 
         ind(i) = array_table(c_index_array1+i-1,target)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, tval1(i), tval2(i))

      enddo

      call three_body_diag(x(isource),x(itarget),
     *                      sval1(1),sval2(1),sval1(2),sval2(2),
     *                      sval1(3),sval2(3),sval1(4),sval2(4), 
     *                      tval1(1),tval2(1),tval1(2),tval2(2),
     *                      tval1(3),tval2(3),epsilon,epsilonb, 
     *                      type(1),type(2),type(3),type(4) ) 

      return
      end

c------------------------------------------------------------------------|
c------------------------------------------------------------------------|

      subroutine three_body_diag(T2,diagel,b1,b2,j1,j2,a1,a2,i1,i2,
     & k1,k2,l1,l2,c1,c2,epsilon,epsilonb,btype,jtype,atype,itype)
      implicit none  
      include 'interpreter.h' 
      include 'parallel_info.h'       


      integer a1,a2,i1,i2,b1,b2,j1,j2,k1,k2,l1,l2,c1,c2,
     & a,b,c,i,j,k,l 
      double precision T2(b1:b2,j1:j2,a1:a2,i1:i2), 
     & diagel(k1:k2,l1:l2,c1:c2) 
      double precision epsilon(*), epsilonb(*)
      integer atype, itype, btype, jtype

      double precision eps, epsa, epsb, epsi, epsj
      double precision val


c-------------Step 1 - contraction---------------------------|

        do a=a1,a2
      
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif

 
c--------------------------------------------|

         do b=b1,b2

          if (btype .eq. mobindex) then
             epsb = epsilonb(b)
          else
             epsb = epsilon(b)
          endif

c--------------------------------------------|

          do i=i1,i2

          if (itype .eq. mobindex) then
            epsi = epsilonb(i)
          else
            epsi = epsilon(i)
          endif

c--------------------------------------------|

           do j=j1,j2
     
          if (jtype .eq. mobindex) then
            epsj = epsilonb(j)
         else
            epsj = epsilon(j)
         endif

c--------------------------------------------|

        eps = epsi + epsj - epsa - epsb

        diagel(j,i,a) = diagel(j,i,a) - T2(b,j,a,i)*T2(b,j,a,i)*eps

 
           enddo
          enddo       
         enddo
        enddo   


c------------ Step 2 - diagonal elements---------------------| 

c           if(me.eq.0) then

             do j=j1,j2
             do i=i1,i2
             do a=a1,a2
c             do b=b1,b2 
c             write(6,*) diagel(j,i,a),'--->',j,i,a
c             enddo
             enddo
             enddo
             enddo 

c             write(6,*) '#########################################'

c           endif  

       return
       end  
