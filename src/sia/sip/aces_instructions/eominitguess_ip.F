c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs diagonalization of small Hbar matrix for the     |
c  construction of the initial guesses for <L| and |R> vectors  |
C                                                               |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine eominitguess_ip(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

c---------------------------------------------------------------|
      implicit none
      include 'interpreter.h'
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
      integer*8 iarray_in, iarray_out, ievec, daddr

      integer i, j, k,n1,n2,m1,m2
      integer array_in,array_out, index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array_out_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
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
c             First array: c_result_array                       |
c---------------------------------------------------------------|
      array_in = op(c_result_array)
C---------------------------------------------------------------|
C             Second array: c_op1_array                         | 
c---------------------------------------------------------------|
      array_out = op(c_op1_array)

C---------------Define type of passed arrays--------------------|
      array_in_type = array_table(c_array_type, array_in)
      array_out_type  = array_table(c_array_type, array_out)
      if (array_in_type .ne. static_array .or.array_out_type 
     &                  .ne. static_array) then
         print *,'Error: Both arrays must be static.'
         call abort_job()
      endif
c---------------------------------------------------------------|
#ifdef ALTIX
      dptr = dshptr
#endif
c---------------Number of indices-------------------------------|
      nindex_in = array_table(c_nindex, array_in)
      nindex_out = array_table(c_nindex,array_out)

c--------------Range of each index------------------------------|

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

      do j = 1, nindex_out
         ind_out(j) = array_table(c_index_array1+j-1,array_out)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
       n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

       m1 = index_table(c_index_size, ind_out(1))  ! pick up length of index
       m2 = index_table(c_index_size, ind_out(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in = get_index_from_base(address_table(array_in), x, 2)
      iarray_out  = get_index_from_base(address_table(array_out), x, 2)

C------Diagonalization intselfs--------------------------------------|

       call initguess_ip(x(iarray_in),1,n1,1,n2,x(iarray_out),
     & 1,m1,1,m2)



       return
       end
C############################################################################
c############################################################################

       subroutine initguess_ip(Hbar,i1,i2,j1,j2,GUESSVEC,a1,a2,b1,b2)
       implicit none 
       include 'daveom.h'
       include 'parallel_info.h'
       integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2,IWRK,LCI,nocc,INFO,
     & ICK,KK,INDX(nroots)
       double precision Hbar(i1:i2,j1:j2),GUESSVEC(a1:a2,b1:b2)
       double precision HTEST(i1:i2,j1:j2),RVEC(i1:i2,j1:j2),
     & LVEC(i1:i2,j1:j2),qqq,fucknorm(nroots),
     & EIG(i2),EIG1(i2),LWRK(8*i2),WR4(i2),WR5(i2),
     & WR1(i2),WR(i2),PPP,testvec(5,5),testvec1(5,5)
c------------------------------
       nocc = i2
       IWRK = 8*nocc
       LCI = nocc

c        if(me.eq.0) write(6,*) i1,i2,'|',j1,j2,'|',a1,a2,'|',b1,b2


c-----------Matrix for diagonalization--------------------------------------|
           DO i=1,LCI
           DO j=1,LCI
           HTEST(i,j) = Hbar(i,j)
           ENDDO
           ENDDO

                        IF(RL.eq.1) then
            call zero(DavRvec_initguess,800*800)
            call zero(DavLvec_initguess,800*800) 
c-----------Diagonalization itselfs-----------------------------------------|
       CALL DGEEV('V','V',LCI,HTEST,LCI,EIG,EIG1,LVEC,LCI,RVEC,LCI,
     & LWRK,IWRK,INFO)


          if(me.eq.0) 
     & write(6,*) 'Initial guesses for EOM-IP eigenvalues'          

          DO i=1,LCI
c          EIG(i) = HTEST(i,i)  
          WR4(i) = EIG(i)
          WR1(i) = EIG(i)

          if(me.eq.0) then
          write(6,*) EIG(i), i  
          endif

          WR5(i) = WR4(i)

          ENDDO
          if(me.eq.0) write(6,*) '|||||||||||||||||||||||||||||||||||||'
         
c--------- Align the guess vectors in the ascending order ------------------|

53    CONTINUE
      ICK=0
      DO KK=1,nocc-1
      IF(EIG(KK+1).lt.EIG(KK)) then
      PPP = EIG(KK)
      EIG(KK) = EIG(KK+1)
      EIG(KK+1) = PPP
      ICK=ICK+1
      ENDIF
      ENDDO
      IF(ICK.NE.0) goto 53

      call zero(INDX,20) 
      DO i=1,nocc
      DO j=1,nocc
      IF(EIG(i).eq.WR1(j)) INDX(i)=j
      ENDDO
      ENDDO

c      if(me.eq.0) write(6,*) EIG(targetroot), INDX(targetroot), '%%%%%%'


c                do i=1,nroots
c          DO j=1,targetroot
c          RVEC(j,i) = 0.d0
c          LVEC(j,i) = 0.d0
c          if(j.eq.i)  then
c          RVEC(j,i) = 1.d0
c          LVEC(j,i)  = 1.d0
c          endif
c          enddo
c                enddo

c------Copy the eigenvalues into the common block----------------------------|
      call zero(DavRoots,800)
      do i=1,nroots
      DavRoots(i) = EIG(targetroot)
      enddo

          do j=1,nroots
          do i=1,nocc
          DavRvec_initguess(i,j) = RVEC(i,INDX(targetroot) )
          DavLvec_initguess(i,j) = LVEC(i,INDX(targetroot) )
          enddo
          enddo
         

c      do i=1,nroots
c      do j=1,nroots
c         qqq = 0.d0
c         do k=1,nocc
c         qqq = qqq + DavRvec_initguess(k,j)*DavLvec_initguess(k,i)
c         enddo
c         if(me.eq.0) write(6,*) qqq,'<--',i,j
c      enddo
c      enddo



          do j=1,nroots
             fucknorm(j) = 0.d0
           do i=1,nocc
           fucknorm(j) = fucknorm(j) + 
     &              DavRvec_initguess(i,j)*DavLvec_initguess(i,j)
           enddo
           enddo


c-------------Renormalization-----------------
          do j=1,nroots
          do i=1,nocc
          DavRvec_initguess(i,j) = 
     &              DavRvec_initguess(i,j) / dsqrt(fucknorm(j))
          DavLvec_initguess(i,j) = 
     &              DavLvec_initguess(i,j) / dsqrt(fucknorm(j))      
          enddo
          enddo
          


          if(me.eq.0) write(6,*) 'Initial guess has been prepared'

                  ENDIF
c-------------Formation of the initial guesses itselfs-----------------------|


c----------------Right eigenvectors------------------------------------------|
          IF(RL.eq.1) then
          do j=1,nroots     
          do i=1,nocc
          GUESSVEC(i,j) = DavRvec_initguess(i,j)/dsqrt(2.d0)
          enddo
          enddo
          if(me.eq.0) write(6,*) 'Copy of right initial guess is done'

          ENDIF
c----------------Left eigenvectors-------------------------------------------|
          IF(RL.eq.0) then
          do j=1,nroots
          do i=1,nocc
          GUESSVEC(i,j) = DavLvec_initguess(i,j)/dsqrt(2.d0)
          enddo
          enddo
          if(me.eq.0) write(6,*) 'Copy of left initial guess is done'
          ENDIF
c----------------------------------------------------------------------------|


       return
       end
