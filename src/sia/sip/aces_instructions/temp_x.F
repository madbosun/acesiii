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
      subroutine temp_x(array_table,
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
      integer*8 iarray, dens_array, indblk, get_block_index
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
c   Get the density array
c----------------------------------------------------------------------------

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. static_array) then
         print *,'Error: density array must be static in ',
     *           'dens2jobarc'
         call abort_job()
      endif

      nindex      = array_table(c_nindex, array)

      if (nindex .ne. 2) then
         print *,'Error: DENS array must have ',
     *           'two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      dens_array = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the type of density
c----------------------------------------------------------------------------
      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         print *,'Error: scalar not sent in dens2jobarc.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in dens2jobarc, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif
      type = INT (scalar_table(sind))

      CALL  DO_THE_WRITING (x(dens_array),inda(1),inda(2),type,
     +                      ibuf,dbuf)

      return
      end
c---------------------------------------------------------------------------

      SUBROUTINE  DO_THE_WRITING (DENS,NDIM1,NDIM2,TYPE,ibuf,dbuf)

      IMPLICIT    NONE

      INCLUDE 'interpreter.h'
      INCLUDE 'int_gen_parms.h'
      INCLUDE 'parallel_info.h'

      INTEGER    NDIM,NDIM1,NDIM2,MU,NU
      INTEGER    ierd,ierd1,ierd2
      INTEGER    TYPE
      INTEGER    ibuf (*)

      CHARACTER*8  JALABEL

      DOUBLE PRECISION  dscale

      DOUBLE PRECISION  DENS (1:NDIM1,1:NDIM2)
      DOUBLE PRECISION  SCR  (1:NDIM1,1:NDIM2)
      DOUBLE PRECISION  SCR1 (1:NDIM1)
      DOUBLE PRECISION  dbuf (*)

      IF (me .EQ. 0) THEN

         IF      (type .EQ.  0) THEN
            JALABEL='HFDENSTY'
         ELSE IF (type .EQ.  1) THEN
            JALABEL='TDENSITY'

         ELSE IF (type .EQ.  2) THEN
            JALABEL='REOMDN01'
         ELSE IF (type .EQ.  3) THEN
            JALABEL='LEOMDN01'
         ELSE IF (type .EQ.  4) THEN
            JALABEL='EXCDEN01'

         ELSE IF (type .EQ.  5) THEN
            JALABEL='REOMDN02'
         ELSE IF (type .EQ.  6) THEN
            JALABEL='LEOMDN02'
         ELSE IF (type .EQ.  7) THEN
            JALABEL='EXCDEN02'

         ELSE IF (type .EQ.  8) THEN
            JALABEL='REOMDN03'
         ELSE IF (type .EQ.  9) THEN
            JALABEL='LEOMDN03'
         ELSE IF (type .EQ. 10) THEN
            JALABEL='EXCDEN03'

         ELSE IF (type .EQ. 11) THEN
            JALABEL='REOMDN04'
         ELSE IF (type .EQ. 12) THEN
            JALABEL='LEOMDN04'
         ELSE IF (type .EQ. 13) THEN
            JALABEL='EXCDEN04'

         ELSE IF (type .EQ. 14) THEN
            JALABEL='REOMDN05'
         ELSE IF (type .EQ. 15) THEN
            JALABEL='LEOMDN05'
         ELSE IF (type .EQ. 16) THEN
            JALABEL='EXCDEN05'

         ELSE IF (type .EQ. 17) THEN
            JALABEL='REOMDN06'
         ELSE IF (type .EQ. 18) THEN
            JALABEL='LEOMDN06'
         ELSE IF (type .EQ. 19) THEN
            JALABEL='EXCDEN06'

         ELSE IF (type .EQ. 20) THEN
            JALABEL='REOMDN07'
         ELSE IF (type .EQ. 21) THEN
            JALABEL='LEOMDN07'
         ELSE IF (type .EQ. 22) THEN
            JALABEL='EXCDEN07'

         ELSE IF (type .EQ. 23) THEN
            JALABEL='REOMDN08'
         ELSE IF (type .EQ. 24) THEN
            JALABEL='LEOMDN08'
         ELSE IF (type .EQ. 25) THEN
            JALABEL='EXCDEN08'

         ELSE IF (type .EQ. 26) THEN
            JALABEL='REOMDN09'
         ELSE IF (type .EQ. 27) THEN
            JALABEL='LEOMDN09'
         ELSE IF (type .EQ. 28) THEN
            JALABEL='EXCDEN09'

         ELSE IF (type .EQ. 29) THEN
            JALABEL='REOMDN10'
         ELSE IF (type .EQ. 30) THEN
            JALABEL='LEOMDN10'
         ELSE IF (type .EQ. 31) THEN
            JALABEL='EXCDEN10'

         ELSE IF (type .EQ. 32) THEN
            JALABEL='REOMDN11'
         ELSE IF (type .EQ. 33) THEN
            JALABEL='LEOMDN11'
         ELSE IF (type .EQ. 34) THEN
            JALABEL='EXCDEN11'

         ELSE IF (type .EQ. 35) THEN
            JALABEL='REOMDN12'
         ELSE IF (type .EQ. 36) THEN
            JALABEL='LEOMDN12'
         ELSE IF (type .EQ. 37) THEN
            JALABEL='EXCDEN12'

         ELSE IF (type .EQ. -1) THEN ! Special for now 5-2011
            JALABEL='AOOVRLAP'
            WRITE (*,*)
            WRITE (*,*) ' Writing atomic overlap integrals to JOBARC '
            WRITE (*,*)

         ELSE

            WRITE (*,*)
            WRITE (*,*) ' No density for this integer label! '
            RETURN

         ENDIF

         IF (intpkg .eq. flocke_package) THEN

             do nu = 1, ndim2
             do mu = 1, ndim1
                scr(mu,nu) = dens(mu,nu)
             enddo
             enddo

          if (type .eq. -1) then 
             do nu = 1, ndim2
             do mu = 1, ndim1
                ierd1 = ibuf (ierdind+mu-1)
                ierd2 = ibuf (ierdind+nu-1)
                dscale = dbuf(iscale_fac+ierd1-1) * 
     *                   dbuf(iscale_fac+ierd2-1)
                dens(ierd1,ierd2) = scr(mu,nu) * dscale
             enddo
             enddo
          else 
             do nu = 1, ndim2
             do mu = 1, ndim1
                ierd1 = ibuf (ierdind+mu-1)
                ierd2 = ibuf (ierdind+nu-1)
                dscale = dbuf(iscale_fac+ierd1-1) * 
     *                   dbuf(iscale_fac+ierd2-1)
                dens(ierd1,ierd2) = scr(mu,nu) / dscale
             enddo
             enddo
          endif 
C
c            do nu = 1, ndim2
c
c               do mu = 1, ndim1
c                  scr1(mu) = dens (mu,nu)
c               enddo
c
c               do mu = 1, ndim1
c                  ierd = ibuf (ierdind+mu-1)
c                  dens (ierd,mu) = scr1(mu) / dbuf (iscale_fac+ierd-1)
c               enddo
c
c            enddo
c
c            do mu = 1, ndim1
c
c               do nu = 1, ndim2
c                  scr1(nu)= dens (mu,nu)
c               enddo
c
c               do nu = 1, ndim2
c                  ierd = ibuf (ierdind+nu-1)
c                  dens (nu,ierd) = scr1(nu) / dbuf (iscale_fac+ierd-1)
c               enddo
c
c            enddo

         ENDIF

         CALL  DESCRAMBLE_MATRIX
     +
     +               (DENS,NDIM1,NFPS,IXSHELLS,NSHELLS)
     +
     +
         NDIM = NDIM1 * NDIM2
         CALL  DPUTREC (1, 'JOBARC', JALABEL, NDIM,DENS)

      ENDIF

      RETURN
      END
