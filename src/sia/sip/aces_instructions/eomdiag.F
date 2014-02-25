c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs a diagonalization of davidson matrix             |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine eomdiag
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
      integer type(mx_array_index) 
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer*8 iarray, ievec, get_index_from_base, daddr
      integer i, j, n1, n2, n3, n4  
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer sind, ind(mx_array_index)
      integer seg, index(mx_array_index), val1(mx_array_index), 
     *        val2(mx_array_index), junk

      integer block, blkndx
      integer find_current_block
      integer*8 get_block_index
      integer stack

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
#ifdef ALTIX
      dptr = dshptr
#endif




c----------------------------------------------------------------------------
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_result_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)

      if (array .eq. 0) then
         print *,'Error: The first argument in return_sval  
     *            must be an array.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index


C--------------Diagonalization itselfs---------------------------------------

      call eom_diag(x(iarray),1,n1,1,n2,iter,nroots,DavRvec,DavLvec,
     & DavRoots)

C----------------------------------------------------------------------------

      return 
      end 
C#############################################################################

      subroutine eom_diag(DAV_Haux,a1,a2,b1,b2,iter,nroots,DavRvec,
     & DavLvec,DavRoots)
      IMPLICIT NONE
      include 'parallel_info.h'

      integer i,j,a1,a2,b1,b2,iter,nroots,IWRK,INFO,LCI,KK,ICK,
     & INDX(800),k
      double precision DAV_Haux(a1:a2,b1:b2),PPP
      double precision EIG(iter*nroots+nroots),EIG1(iter*nroots+nroots),
     & VECTOR(iter*nroots+nroots,iter*nroots+nroots),
     & VECTOR1(iter*nroots+nroots,iter*nroots+nroots),
     & LWRK(8*(iter*nroots+nroots)),WR5(iter*nroots+nroots),
     & HTEST(iter*nroots+nroots,iter*nroots+nroots),
     & TV(iter*nroots+nroots),WR4(iter*nroots+nroots),
     & DavRvec(800,800),DavLvec(800,800),DavRoots(800)
      
      IWRK = 8*(iter*nroots+nroots) 
      LCI = iter*nroots+nroots

c-----------Matrix for diagonalization----------------------------------------

           DO i=1,LCI
           DO j=1,LCI
           HTEST(i,j) = DAV_Haux(i,j)
c           if(me.eq.0) write(6,*) HTEST(i,j),'!!!!!!!!!!!!!!!!!',i,j 
           ENDDO
           ENDDO
c-----------Diagonalization itselfs-------------------------------------------
       CALL DGEEV('V','V',LCI,HTEST,LCI,EIG,EIG1,VECTOR1,LCI,VECTOR,LCI,
     & LWRK,IWRK,INFO)

          if(me.eq.0) write(6,*) 'On iteration',iter,'eigenvalues are:' 
          DO i=1,LCI
          WR4(i) = EIG(i)
          WR5(i) = WR4(i)
          if(me.eq.0) write(6,*) EIG(i),EIG1(i),'~~~~~~~>',i
          ENDDO
c--------- Align the eigenvalues in the discending order------------------|

63    CONTINUE

      ICK=0
      DO KK=1,LCI-1
      IF( (WR4(KK+1)).lt.(WR4(KK)) ) then
      PPP = WR4(KK)
      WR4(KK) = WR4(KK+1)
      WR4(KK+1) = PPP
      ICK=ICK+1
      ENDIF
      ENDDO
      IF(ICK.NE.0) goto 63

      call zero(INDX,800)
      DO j=1,LCI
      DO i=1,LCI
      IF(WR4(i).eq.EIG(j)) INDX(i)=j
      ENDDO
      ENDDO

c-------Copy roots of small matrix to the common block------------------------
      call zero(DavRoots,800)

      do i=1,LCI
      DavRoots(i) = WR4(i)
      if(me.eq.0) write(6,*) DavRoots(i), '/\/\/\/\/\/\/\/\/',i  
      enddo

c-------Copy left and right auxiliary vectors to the common block-------------

      call zero(DavRvec,800*800)
      call zero(DavLvec,800*800)

         DO i=1,LCI
         DO j=1,LCI
         DavRvec(i,j) = VECTOR(i, INDX(j) )
         DavLvec(i,j) = VECTOR1(i,INDX(j) ) 
         ENDDO
         ENDDO

c      do i=1,LCI
c      do j=1,LCI
c         ppp = 0.d0
c         do k=1,LCI
c         ppp = ppp+DavLvec(k,i)*DavRvec(k,j)
c        enddo
c         if(me.eq.0) write(6,*) ppp,'<--',i,j
c      enddo
c      enddo



      return
      end
c-----------------------------------------------------------------------------
