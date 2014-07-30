c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It does create overlap table of new eom roots and roots from |
c             previous iteration                                |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine eomovertable
     +                 (array_table,
     *                  narray_table, index_table,
     *                  nindex_table, segment_table, nsegment_table,
     *                  block_map_table, nblock_map_table,
     *                  scalar_table, nscalar_table,
     *                  address_table, op)

c--------------------------------------------------------------------------
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'epsilon.h'
      include 'daveom.h'
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
      integer*8 indblk1, indblk2, get_block_index,iarray,daddr
      integer stack
      integer nsend
      integer array_type, evec_type
      integer comm
      integer sind, ind(mx_array_index)
      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2,n1,n2
      integer*8 addr, get_index_from_base
      double precision x(1)


      double precision overmax_index

#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------|
c   Locate the data for the input matrix which is to be analyses.           |
c   First array: set of overlaps for current root                           |
c---------------------------------------------------------------------------|

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
c---------------------------------------------------------------------------|
c               Selection of maximal root                                   |
c---------------------------------------------------------------------------|

        call max_overlap2( x(iarray), n1, n2 )

c---------------------------------------------------------------------------|
c---------------------------------------------------------------------------|



      return
      end
c===========================================================================|
c===========================================================================|
      subroutine max_overlap2(overlap,i1,i2)
      implicit none
      include 'parallel_info.h'
      include 'daveom.h'
      integer i1,i2,i,j,kk,overmax,izz,iyy
      double precision overlap(i1,i2),overmax_index,qq,qq1,
     & root_tmp(400)
c-----------------------------------------------------------|
c This subroutine does search of maximal overlap            |
c-----------------------------------------------------------|

            izz = nroots*iter+nroots  

            DO i=1,nroots              

            overmax = 1
            DO kk=1,izz
            qq = overlap(overmax,i)
            qq1 = overlap(kk,i)

c            IF(RL.eq.0) then                   
c            if(me.eq.0) write(6,*) dabs(qq1),'****',i,'****',kk
c            endif                 

                  IF( dabs(qq1) .gt. dabs(qq) ) then
                  overmax = kk
                  ENDIF

            ENDDO

                 IF(RL.eq.1) then
                 DavOverRmax(i) = overmax
                 if(me.eq.0) write(6,*) DavOverRmax(i),'***===***',i
                 ENDIF

                 IF(RL.eq.0) then
                 DavOverLmax(i) = overmax 
                 if(me.eq.0) write(6,*) DavOverLmax(i),'***$$$***',i 
                 ENDIF  

            ENDDO 

c--------Reordering of the eigenvalues according to maximal overlap----
              IF(RL.eq.1) then
            DO i=1,izz  
            root_tmp(i) = DavRoots(i)
            ENDDO


            DO i=1,izz 
            overmax = DavOverRmax(i)
            DavRoots(i) = root_tmp(overmax)
            ENDDO
              ENDIF 
c----------------------------------------------------------------------
                 
           


c            if(me.eq.0) write(6,*) overmax,'#####',overlap(overmax)

              

      return
      end 
C---------------------------------------------------------------------------|





