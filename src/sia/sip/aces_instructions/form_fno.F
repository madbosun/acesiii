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
      subroutine form_fno(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Both Array1 and Array2 must be "static" arrays.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'trace.h'
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

      integer*8 iarray, ievec, get_index_from_base
      integer i, j, n
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer ind(mx_array_index)
      integer junk

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

c----------------------------------------------------------------------------
c   Locate the data for the input matrix.
c---------------------------------------------------------------------------

      array = op(c_result_array)
      evec_array = op(c_op1_array)
      if (array .eq. 0 .or. evec_array .eq. 0) then
         print *,'Error: EIG routine requires 2 array arguments.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      array_type = array_table(c_array_type, array)
      evec_type  = array_table(c_array_type, evec_array)
      if (array_type .ne. static_array .or. 
     *    evec_type .ne. static_array) then
         print *,'Error: Both arrays in EIG routine must be static.'
         call abort_job()
      endif

#ifdef ALTIX
      dptr = dshptr
#endif

      nindex = array_table(c_nindex, array)
      nindex_evec = array_table(c_nindex, evec_array)
      if (nindex .ne. nindex_evec) then
         print *,'Error: Both arrays in EIG must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex,' indices.'
         print *,'Second array has ',nindex_evec,' indices.'
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray = get_index_from_base(address_table(array), x, 2)
      ievec  = get_index_from_base(address_table(evec_array), x, 2)

      n = index_table(c_index_size, ind(1))  ! pick up length of index
      call do_formfno(x(iarray), x(ievec), n, ind(1), ind(2))

      return
      end

      subroutine do_formfno(SPab, Uab, n, n1, n2)
      include 'int_gen_parms.h'
      parameter (max_basis = 2000)
      integer i, j, n, n1, n2, k
      double precision SPab(n,n)
      double precision Uab(n,n)  
      double precision utemp(n,n)  
      double precision thresh, maxp, minp, temp, temp1, dtmp
      integer ndrop
      integer ndropmo_a
      integer modrop_a(max_basis)
      logical harddrop

      open(66,file='summary.out') 
      if (me .eq. 0) write(66,5)  

      thresh = 0.00005d0
      harddrop = .false.
      ndrop = 0 
      minp = SPab(1,1)
      maxp = SPab(n,n)

      if (me .eq. 0) write(66,40) minp,maxp

      call igetrec(1, 'JOBARC', 'OCCUPYA0', 1, nalpha_occupied)
      CALL IGETREC(20,'JOBARC','NUMDROPA',1,NDROPMO_A)
      if (ndropmo_a .gt. 0) then
          harddrop = .true.
          CALL IGETREC(20,'JOBARC','MODROPA',NDROPMO_A,MODROP_A)
          do i = 1,ndropmo_a
              if(modrop_a(i) .gt. nalpha_occupied) then
                  ndrop = ndrop + 1
              end if
          end do
      else
          harddrop = .false.
          ndrop = 0
          if (me .eq. 0) write(66,21) thresh 
      end if
C      if(me.eq.0) write(66,*) fno_keep

      do j = 1,n
         if (harddrop .and. j.le.ndrop) then
C            if (me .eq. 0) write(66,11) j+nalpha_occupied, SPab(j,j),'D'
      if (me .eq. 0) write(66,11) n+nalpha_occupied-j+1, SPab(j,j),'D'
            do i = 1, n
               Uab(i,j) = 0.0 
            enddo
         else if (.not.harddrop .and. dabs(SPab(j,j)) .lt. thresh) then 
            ndrop = ndrop + 1 
C            if (me .eq. 0) write(66,11) j+nalpha_occupied, SPab(j,j),'D'
      if (me .eq. 0) write(66,11) n+nalpha_occupied-j+1, SPab(j,j),'D'
            do i = 1, n
               Uab(i,j) = 0.0 
            enddo
         else
C            if (me .eq. 0) write(66,11) j+nalpha_occupied, SPab(j,j),''
      if (me .eq. 0) write(66,11) n+nalpha_occupied-j+1, SPab(j,j),' '
         endif 
      enddo

      if (me .eq. 0) write(66,30) ndrop, n  
      if (me .eq. 0) 
     *   write(66,31) n+nalpha_occupied-ndrop,n+nalpha_occupied

5     format(5x,'DETERMINING MBPT(2) FROZEN NATURAL ORBITALS',/)  
40    format(/5x,'MIN/MAX virtual occupancy=',F12.8,'/',F12.8)  
20    format(5x,'DETERMINING FNO SUBSPACE BASED ON % DROP: ',I3,'%')  
21    format(5x,'DETERMINING FNO SUBSPACE BASED ON THRESHOLD: ',F12.8)  
999   format(I4,100F10.7)
11    format(5x,I3,1F14.8,2x,A)
30    format(5x,'Dropping ',I8,' natural orbitals out of ',I8,
     *       ' total virtuals')
31    format(5x,'If this is the first run, use DROPMO=',I8,'-',I8)

      return
      end
