C  Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
C 
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C  USA
      subroutine form_fock_rohf(array_table,
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

      integer array_type, nindex, inda(2), indb(2)    

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 iarray, farray, indblk, get_block_index
      integer*8 ibarray, fbarray 
      integer*8 integral_scr
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack

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

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c----------------------------------------------------------------------------
c   Get the alpha Fock array.
c----------------------------------------------------------------------------

      array = 9 ! Fock_a for a placeholder   
      array_type = array_table(c_array_type, array)

      if (array_type .ne. static_array) then
         print *,'Error: Fock arrays in integral calc
     *            must be static.'
         call abort_job()
      endif

      nindex      = array_table(c_nindex, array)

      if (nindex .ne. 2) then
         print *,'Error: Fock array must have ',
     *           'two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         indb(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      farray = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the beta Fock array.
c----------------------------------------------------------------------------

      array = 10 ! Fock_b for a placeholder   
      array_type = array_table(c_array_type, array)

      if (array_type .ne. static_array) then
         print *,'Error: Fock arrays in integral calc
     *            must be static.'
         call abort_job()
      endif

      nindex      = array_table(c_nindex, array)

      if (nindex .ne. 2) then
         print *,'Error: Fock array must have ',
     *           'two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         indb(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      fbarray = get_index_from_base(address_table(array), x, 2)

      if (.not. dryrun) 
     *   call fock_rohf(1, indb(1), 1,indb(2), x(farray),x(fbarray))

      return
      end

      subroutine fock_rohf(s1,nend1,s2,nend2,Fapb,Famb) 
      implicit none
      include 'int_gen_parms.h'
      include 'parallel_info.h'

      integer s1, s2, nend1, nend2  
      integer i, j, k, l 
      integer nsingle, ndouble 
      double precision Fapb(1:nend1,1:nend2)   
      double precision Famb(1:nend1,1:nend2)   
      double precision Fa(1:nend1,1:nend2)   

      if (nalpha_occupied .eq. nbeta_occupied) then 
         ndouble = nalpha_occupied 
         nsingle = 0 
      else 
        ndouble = min(nalpha_occupied,nbeta_occupied) 
        nsingle = abs(nalpha_occupied-nbeta_occupied) 
      endif 
c
c Zero out the Fock_a(ROHF) array 
c
      DO I = 1, nbasis 
      DO J = 1, nbasis 
         FA(I,J) = 0.0d0  
      ENDDO 
      ENDDO 
c
c DOCC-DOCC block (arbitrary and it is set to -1/2FA + 3/2FB)
c
      Do I = 1, ndouble   
         Do J = 1, ndouble  
            FA(I,J) = (FAPB(I,J) + FAMB(I,J))/2.0d0 
         Enddo
      Enddo
c
cSOCC-SOCC block (arbitrary and it is set to 1/2FA + 1/2FB)
c
      Do I = ndouble+1, ndouble+nsingle  
         Do J = ndouble+1, ndouble+nsingle  
            FA(I,J) = (FAPB(I,J) + FAMB(I,J))/2.0d0 
         Enddo
      Enddo
c
cVIRT-VIRT block (arbitrary and it is set to 3/2FA - 1/2FB)
c
      Do I = ndouble+nsingle+1, nbasis  
         Do J = ndouble+nsingle+1, nbasis  
            FA(I,J) = (FAPB(I,J) + FAMB(I,J))/2.0d0 
         Enddo
      Enddo
c
cDOCC-SOCC block (arbitrary and it is set to FB)
c
      DO I = 1, ndouble  
         Do J = ndouble+1, ndouble+nsingle  
            FA(I,J) = FAMB(I,J)
         Enddo
      Enddo
c
      DO J = 1, ndouble  
         Do I = ndouble+1, ndouble+nsingle  
            FA(I,J) = FAMB(I,J)
         Enddo
      Enddo
c
cVIRT-DOCC block (arbitrary and it is set to 1/2(FB+FA))
c
      DO I = ndouble+nsingle+1, nbasis  
         Do J = 1, ndouble  
            FA(I,J) = (FAPB(I,J) + FAMB(I,J))/2.0d0 
         Enddo
      Enddo
c
      DO J = ndouble+nsingle+1, nbasis  
         Do I = 1, ndouble  
            FA(I,J) = (FAPB(I,J) + FAMB(I,J))/2.0d0 
         Enddo
      Enddo
c
cVIRT-SOCC block (arbitrary and it is set to FA)
c
      DO I = ndouble+nsingle+1, nbasis  
         Do J = ndouble+1, ndouble+nsingle  
            FA(I,J) = FAPB(I,J) 
         Enddo
      Enddo
c
      DO J = ndouble+nsingle+1, nbasis  
         Do I = ndouble+1, ndouble+nsingle  
            FA(I,J) = FAPB(I,J) 
         Enddo
      Enddo

100   continue 

      DO I = 1, nbasis 
      DO J = 1, nbasis 
         if (dabs(FA(I,J)) .GT. 100d0) 
     *       write(6,*) ' LARGE F ELEMENT', i,j,fa(i,j)
         FAMB(I,J) = FA(I,J) 
         FAPB(I,J) = FA(I,J) 
      ENDDO 
      ENDDO 
c     Do I = 1, nbasis    
c     Do J = 1, nbasis 
c        write(6,*) 'F', i, j, FAPB(I,J) , FAMB(I,J)  
c     enddo 
c     enddo 

      return 
      end 
