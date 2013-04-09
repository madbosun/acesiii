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
      subroutine write_coefficients(array_table,
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
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, array 
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table), get_index_from_base
      integer*8 get_block_data_index
      integer find_current_block, block 

      integer array_type, nindex, inda(2), indb(2), indc(2),indd(2)  

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 iarray, ca_array, cb_array, indblk, get_block_index
      integer*8 fa_array, fb_array
      integer*8 ibarray, fbarray 
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack
      integer sind,type  

      double precision x(1)
      double precision dbuf(1), temp 
      integer ibuf(1)
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

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c----------------------------------------------------------------------------
c   Get the CA
c----------------------------------------------------------------------------

      array = 2
      array_type = array_table(c_array_type, array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo
      ca_array = get_index_from_base(address_table(array), x, 2)

      array = 3
      array_type = array_table(c_array_type, array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         nb = array_table(c_index_array1+i-1,array)
         indb(i) = index_table(c_index_size, nb)  ! pick up length of index
      enddo
      cb_array = get_index_from_base(address_table(array), x, 2)

      array = 9
      array_type = array_table(c_array_type, array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         nc = array_table(c_index_array1+i-1,array)
         indc(i) = index_table(c_index_size, nc)  ! pick up length of index
      enddo
      fb_array = get_index_from_base(address_table(array), x, 2)

      array = 10
      array_type = array_table(c_array_type, array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         nd = array_table(c_index_array1+i-1,array)
         indd(i) = index_table(c_index_size, nd)  ! pick up length of index
      enddo
      fa_array = get_index_from_base(address_table(array), x, 2)


      CALL  WRITE_C (x(ca_array),inda(1),inda(2),
     +               x(cb_array),indb(1),indb(2),
     +               x(fa_array),indc(1),indc(2),
     +               x(fb_array),indd(1),indd(2))

      return
      end
c---------------------------------------------------------------------------

      SUBROUTINE  WRITE_C (CA,ADIM1,ADIM2,
     +                     CB,BDIM1,BDIM2,
     +                     FA,CDIM1,CDIM2,
     +                     FB,DDIM1,DDIM2 )

      IMPLICIT    NONE

      INCLUDE 'interpreter.h'
      INCLUDE 'int_gen_parms.h'
      INCLUDE 'parallel_info.h'

      INTEGER    ADIM1,ADIM2,BDIM1,BDIM2
      INTEGER    CDIM1,CDIM2,DDIM1,DDIM2
      INTEGER    TYPE,P,Q,MU

      DOUBLE PRECISION  CA (1:ADIM1,1:ADIM2)
      DOUBLE PRECISION  CB (1:BDIM1,1:BDIM2)
      DOUBLE PRECISION  FA (1:CDIM1,1:CDIM2)
      DOUBLE PRECISION  FB (1:DDIM1,1:DDIM2)

      IF (me .EQ. 0) THEN

         OPEN(UNIT=500,FILE='COEFF',STATUS='UNKNOWN')
         REWIND (500)

         DO P = 1, ADIM1
         DO Q = 1, ADIM2
            WRITE(500,*) CA (P,Q)
         ENDDO
         ENDDO

         DO P = 1, BDIM1
         DO Q = 1, BDIM2
            WRITE(500,*) CB (P,Q)
         ENDDO
         ENDDO
 
         DO P = 1, CDIM1
         DO Q = 1, CDIM2
            WRITE(500,*) FA (P,Q)
         ENDDO
         ENDDO

         DO P = 1, DDIM1
         DO Q = 1, DDIM2
            WRITE(500,*) FB (P,Q)
         ENDDO
         ENDDO

         CLOSE(500,STATUS='KEEP')

      ENDIF

      RETURN
      END
