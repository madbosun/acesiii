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
      subroutine prestrict

     +                     (array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

c--------------------------------------------------------------------------
c   
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'pst_functions.h'
      include 'int_gen_parms.h'
      include 'restrict_pardo.h'
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
      integer evec_array1, evec_type1, nindex_evec1
      integer sind, sind1, ind(mx_array_index)
      integer seg, index(mx_array_index), val1(mx_array_index),
     *        val2(mx_array_index), junk

      integer block, blkndx
      integer find_current_block
      integer*8 get_block_index
      integer stack

      integer pcount

      integer nakuldebug

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
c   Locate the scalars. --> bproc and eproc  
c---------------------------------------------------------------------------

      evec_array  = op(c_result_array)
      evec_type   = array_table(c_array_type, evec_array)
      nindex_evec = array_table(c_nindex, evec_array)

      evec_array1 = op(c_op1_array)
      evec_type1  = array_table(c_array_type, evec_array1)
      nindex_evec1= array_table(c_nindex, evec_array1)


      if (evec_type .ne. scalar_value) return
      if (evec_type1.ne. scalar_value) return

      sind =  scalar_table(array_table(c_scalar_index, evec_array))
      sind1=  scalar_table(array_table(c_scalar_index, evec_array1))

c     if (evec_array .lt. 1 .or. evec_array .gt. nscalar_table) then
c        print *,'Scalar table index 1 out of range in prestrict, ',
c    *           'line ',current_line
c        print *,'Index for array ',evec_array,'is',sind,' should be ',
c    *           'between 1 and ',nscalar_table
c        call abort_job()
c     endif

c     if (evec_array1.lt. 1 .or. evec_array1 .gt. nscalar_table) then
c        print *,'Scalar table index 2 out of range in prestrict, ',
c    *           'line ',current_line
c        print *,'Index for array ',evec_array,'is',sind1,' should be ',
c    *           'between 1 and ',nscalar_table
c        call abort_job()
c     endif

c     write(6,*) ' SIND :', sind, sind1 

      pcount = 0
      do i = 1, nprocs
         if (pst_get_company(i-1) .eq. io_company_id) then
c           write(6,*) ' IO SERVER :', i
            do_prestrict(i) = 0
         else
            pcount = pcount + 1
            if (pcount .ge. sind .and. pcount .le. sind1) then
               do_prestrict(i) = 1
            else
               do_prestrict(i) = 0
            endif
c           write(6,*) ' WORKER :', i
         endif
c         print *, i, do_prestrict(i), ' in prestrict'
      enddo

c       nakuldebug = 0

c       do i = 1, nprocs
c           if (do_prestrict(i) .eq. 1) 
c     *          nakuldebug = nakuldebug +1
c       enddo

c      print *, me, ' ', nakuldebug
c       if (nakuldebug .ne. 24 .and. 
c     *      nakuldebug .ne. 4. and.
c     *      nakuldebug .ne. 17)
c     *      print *, ' corrupted do_prestrict table'



c     write(6,*) ' DO_PRESTRICT:', (do_prestrict(i), i=1, nprocs) 


      return
      end


