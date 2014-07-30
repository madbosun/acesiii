c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs construction of new approximations to |R> and    |
c  and <L| vectors using Davodson's vectors and eigenvectors    |
c  of the small matrix                                          |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine eomnewvec(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

c---------------------------------------------------------------|
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
      include 'daveom.h' 
#ifdef ALTIX
      include 'sheap.h'
#endif
c---------------------------------------------------------------|

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
      integer*8 indblk1, indblk2, get_block_index
      integer stack
      integer nsend

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer*8 addr, get_index_from_base
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------|

c---------------------------------------------------------------|
c   Detrmine the parameters of the first array: c_result_array  |
c---------------------------------------------------------------|

      array = op(c_result_array)
      nindex = array_table(c_nindex, array)

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         sindex(i) = index
         call get_index_segment(index, seg, segment_table,
     *                           nsegment_table, index_table,
     *                           nindex_table, fop1(i), fop2(i))
      enddo


      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk1 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk1 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif
c---------------------------------------------------------------|
c   Detrmine the parameters of the second array: c_op1_array    |
c---------------------------------------------------------------|
      array = op(c_op1_array)
      nindex = array_table(c_nindex, array)

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         findex(i) = index
         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, sop1(i), sop2(i))
      enddo


      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk2 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk2 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif
C---------------------------------------------------------------|
           
C--------------Check type of contraction------------------------|

c-----Contraction with 1 index-Davidson vectors contribution ---|
         if (chunk_type.eq.1) then 
         call make_eom_newvec_1( x(indblk1),fop1(1),fop2(1),
     & fop1(2),fop2(2),fop1(3),fop2(3),x(indblk2),
     & sop1(1),sop2(1),sop1(2),sop2(2) ) 
         endif 
c-----Contraction with 1 index-correction vectors contribution --|
         if (chunk_type.eq.2) then                   
         call make_eom_newvec_2( x(indblk1),fop1(1),fop2(1),
     & fop1(2),fop2(2),x(indblk2),sop1(1),sop2(1),sop1(2),sop2(2) )   
         endif
c-----Contraction with 3 indices- Davidson vectors contribution--|
         if (chunk_type.eq.3) then
       call make_eom_newvec_3( x(indblk1),fop1(1),fop2(1),
     & fop1(2),fop2(2),fop1(3),fop2(3),fop1(4),fop2(4),fop1(5),fop2(5),
     & x(indblk2),sop1(1),sop2(1),sop1(2),sop2(2),sop1(3),sop2(3),
     & sop1(4),sop2(4) ) 
         endif     
c-----Contraction with 3 indices- correction vectors contribution--|
         if (chunk_type.eq.4) then
       call make_eom_newvec_4( x(indblk1),fop1(1),fop2(1),
     & fop1(2),fop2(2),fop1(3),fop2(3),fop1(4),fop2(4),
     & x(indblk2),sop1(1),sop2(1),sop1(2),sop2(2),sop1(3),sop2(3),
     & sop1(4),sop2(4) )
         endif



       return
       end
C############################################################################
C############################################################################
c        One-index component                                                |
c---------------------------------------------------------------------------|
       subroutine make_eom_newvec_1(DAV_VEC,a1,a2,b1,b2,c1,c2,
     & DAV_NEWVEC,d1,d2,f1,f2 ) 
       implicit none
       include 'daveom.h'   
       integer a1,a2,b1,b2,c1,c2,d1,d2,f1,f2,i,j,a,kk,ll,zz,qq,ist,
     & izz
       double precision DAV_VEC(a1:a2,b1:b2,c1:c2),
     & DAV_NEWVEC(d1:d2,f1:f2)

       ist = nroots*iter
       izz = nroots*iter + nroots  
c--------Davidson's vectors contribution------------------------------------|
       
cccc           DO zz=1,izz !loop over roots of output array
            DO i=a1,a2 !loop over hole molecular index
             qq = 0
             DO ll=c1,c2 !loop over iterations
              DO kk=b1,b2 !loop over roots of input array
              qq = qq + 1
                 IF(RL.eq.1) then
              DAV_NEWVEC(i,aux_vec_num) = DAV_NEWVEC(i,zz) + 
     & DAV_VEC(i,kk,ll)*DavRvec(qq,aux_vec_num) 
                 ELSE
              DAV_NEWVEC(i,aux_vec_num) = DAV_NEWVEC(i,zz) +
     & DAV_VEC(i,kk,ll)*DavLvec(qq,aux_vec_num)
                 ENDIF
              ENDDO
             ENDDO 
            ENDDO
cccc           ENDDO
      
c---------------------------------------------------------------------------|
      return
      end   
C############################################################################
C############################################################################
       subroutine make_eom_newvec_2(DAV_VEC,a1,a2,b1,b2,DAV_NEWVEC,
     & d1,d2,f1,f2 )
       implicit none
       include 'daveom.h'
       integer a1,a2,b1,b2,c1,c2,d1,d2,f1,f2,i,j,a,kk,ll,zz,qq,ist,
     & izz
       double precision DAV_VEC(a1:a2,b1:b2),DAV_NEWVEC(d1:d2,f1:f2)

       ist = nroots*iter
       izz = nroots*iter + nroots
c--------Correction vectors contribution------------------------------------|
    
cccc           DO zz=1,izz !loop over roots of output array
            DO i=a1,a2 !loop over hole molecular index
             DO kk=b1,b2 !loop over roots of input array
                 IF(RL.eq.1) then
              DAV_NEWVEC(i,aux_vec_num) = DAV_NEWVEC(i,zz) +
     & DAV_VEC(i,kk)*DavRvec(kk+ist,aux_vec_num)
                 ELSE
              DAV_NEWVEC(i,aux_vec_num) = DAV_NEWVEC(i,zz) +
     & DAV_VEC(i,kk)*DavLvec(kk+ist,aux_vec_num)
                 ENDIF
             ENDDO
            ENDDO
cccc           ENDDO
       
c---------------------------------------------------------------------------|
      return
      end
C############################################################################
C############################################################################
c                Three-index component                                      |
c---------------------------------------------------------------------------|
       subroutine make_eom_newvec_3(DAV_VEC,a1,a2,b1,b2,c1,c2,d1,d2,
     & f1,f2,DAV_NEWVEC,g1,g2,h1,h2,k1,k2,l1,l2)
       implicit none
       include 'daveom.h'
       integer a1,a2,b1,b2,c1,c2,d1,d2,f1,f2,i,j,a,kk,ll,l1,l2,h1,h2,g1,
     & g2,k1,k2,zz,qq,izz
       double precision DAV_VEC(a1:a2,b1:b2,c1:c2,d1:d2,f1:f2),
     & DAV_NEWVEC(g1:g2,h1:h2,k1:k2,l1:l2)

       izz = nroots*iter + nroots
c--------Davidson's vectors contribution------------------------------------|

cccc           DO zz=1,izz !loop over roots of output array
            DO i=a1,a2 !loop over hole molecular index
            DO j=b1,b2 !loop over hole molecular index
            DO a=c1,c2 !loop over particle molecular index
             qq = 0
             DO ll=f1,f2 !loop over iterations
              DO kk=d1,d2 !loop over roots of input array
              qq = qq + 1
                 IF(RL.eq.1) then
              DAV_NEWVEC(i,j,a,aux_vec_num) = DAV_NEWVEC(i,j,a,zz) +
     & DAV_VEC(i,j,a,kk,ll)*DavRvec(qq,aux_vec_num)
                 ELSE
              DAV_NEWVEC(i,j,a,aux_vec_num) = DAV_NEWVEC(i,j,a,zz) +
     & DAV_VEC(i,j,a,kk,ll)*DavLvec(qq,aux_vec_num)
                 ENDIF
              ENDDO
             ENDDO
            ENDDO
            ENDDO
            ENDDO
cccc           ENDDO
c---------------------------------------------------------------------------|
      return
      end
C############################################################################
C############################################################################
c                Three-index component                                      |
c---------------------------------------------------------------------------|
       subroutine make_eom_newvec_4(DAV_VEC,a1,a2,b1,b2,c1,c2,d1,d2,
     & DAV_NEWVEC,g1,g2,h1,h2,k1,k2,l1,l2)
       implicit none
       include 'daveom.h'
       integer a1,a2,b1,b2,c1,c2,d1,d2,f1,f2,i,j,a,kk,ll,l1,l2,h1,h2,g1,
     & g2,k1,k2,ist,zz,izz
       double precision DAV_VEC(a1:a2,b1:b2,c1:c2,d1:d2),
     & DAV_NEWVEC(g1:g2,h1:h2,k1:k2,l1:l2)


       ist = nroots*iter
       izz = nroots*iter + nroots
c--------Correction vectors contribution--------------------------------------|

cccc           DO zz=1,izz !loop over roots of output array
            DO i=a1,a2 !loop over hole molecular index
            DO j=b1,b2 !loop over hole molecular index
            DO a=c1,c2 !loop over particle molecular index
              DO kk=d1,d2 !loop over roots of input array
                  IF(RL.eq.1) then
              DAV_NEWVEC(i,j,a,aux_vec_num) = DAV_NEWVEC(i,j,a,zz) +
     & DAV_VEC(i,j,a,kk)*DavRvec(kk+ist,aux_vec_num)
                  ELSE
              DAV_NEWVEC(i,j,a,aux_vec_num) = DAV_NEWVEC(i,j,a,zz) +
     & DAV_VEC(i,j,a,kk)*DavLvec(kk+ist,aux_vec_num)
                  ENDIF
              ENDDO
            ENDDO
            ENDDO
            ENDDO
cccc           ENDDO
c---------------------------------------------------------------------------|
      return
      end
c----------------------------------------------------------------------------





