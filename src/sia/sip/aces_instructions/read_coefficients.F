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
      subroutine  read_coefficients(array_table,
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
      include 'epsilon.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer read_type
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

      integer array_type, nindex, inda(2), indb(2), indc(2), indd(2)
      integer inde(2), indf(2), ne, nf, nindex_evec

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk, evec_array, evec_type
      integer*8 iarray, ca_array, cb_array, indblk, get_block_index
      integer*8 fa_array, fb_array, ea_array, eb_array 
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
      fa_array = get_index_from_base(address_table(array), x, 2)

      array = 10
      array_type = array_table(c_array_type, array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         nd = array_table(c_index_array1+i-1,array)
         indd(i) = index_table(c_index_size, nd)  ! pick up length of index
      enddo
      fb_array = get_index_from_base(address_table(array), x, 2)

      evec_array  = op(c_result_array)
      evec_type   = array_table(c_array_type, evec_array)
      nindex_evec = array_table(c_nindex, evec_array)

      if (evec_type .ne. scalar_value) return

      if (evec_array .lt. 1 .or. evec_array .gt. narray_table) then
         print *,'Error: Invalid scalar in return_sval, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif
      
      sind =  array_table(c_scalar_index, evec_array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in return_sval, ',
     *           'line ',current_line
         print *,'Index for array ',evec_array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      if ( (abs (scalar_table(sind) - 1.0)) .lt. 0.0001 ) read_type = 1
      if ( (abs (scalar_table(sind) - 2.0)) .lt. 0.0001 ) read_type = 2

      CALL   READ_C (x(ca_array),inda(1),inda(2),
     +               x(cb_array),indb(1),indb(2),
     +               x(fa_array),indc(1),indc(2),
     +               x(fb_array),indd(1),indd(2),
     +               epsilon,
     +               epsilonb,read_type)

      return
      end
c---------------------------------------------------------------------------

      SUBROUTINE   READ_C (CA,ADIM1,ADIM2,
     +                     CB,BDIM1,BDIM2,
     +                     FA,CDIM1,CDIM2,
     +                     FB,DDIM1,DDIM2,
     +                     EA,
     +                     EB,read_type)

      IMPLICIT    NONE

      INCLUDE 'interpreter.h'
      INCLUDE 'int_gen_parms.h'
      INCLUDE 'parallel_info.h'

      INTEGER    read_type
      INTEGER    ADIM1,ADIM2,BDIM1,BDIM2
      INTEGER    CDIM1,CDIM2,DDIM1,DDIM2
      INTEGER    EDIM1,FDIM1,NALPHA,NBETA
      INTEGER    ALPHA_ORB,BETA_ORB
      INTEGER    TYPE,P,Q,MU

      LOGICAL    ASWAP,BSWAP,BREAK_SYM

      DOUBLE PRECISION  DROOT2

      DOUBLE PRECISION  CA (1:ADIM1,1:ADIM2)
      DOUBLE PRECISION  CB (1:BDIM1,1:BDIM2)
      DOUBLE PRECISION  FA (1:CDIM1,1:CDIM2)
      DOUBLE PRECISION  FB (1:DDIM1,1:DDIM2)
      DOUBLE PRECISION  EA (1:CDIM1)
      DOUBLE PRECISION  EB (1:DDIM1)

      DOUBLE PRECISION  TMP (1:BDIM1,1:BDIM2)
C
C
C   Reading the orbitals to perform the exact N-1 calculation
C
C
      IF (read_type .EQ. 1) THEN
      IF (me .EQ. 0) THEN
C
C
C   Read the coefficients for the N electron system!
C
C
         OPEN(UNIT=500,FILE='COEFF',STATUS='UNKNOWN')
         REWIND (500)

         DO MU = 1, ADIM1
         DO Q = 1, ADIM2
            READ (500,*) CA (MU,Q)
         ENDDO
         ENDDO

         DO MU = 1, BDIM1
         DO Q = 1, BDIM2
            READ (500,*) CB (MU,Q)
         ENDDO
         ENDDO

         DO P = 1, CDIM1
         DO Q = 1, CDIM2
            READ (500,*) FA (P,Q)
         ENDDO
         ENDDO

         DO P = 1, DDIM1
         DO Q = 1, DDIM2
            READ (500,*) FB (P,Q)
         ENDDO
         ENDDO

         DO P = 1, CDIM1
            EA (P) = FA (P,P)
         ENDDO 

         DO Q = 1, DDIM1
            EB (Q) = FB (Q,Q)
         ENDDO
C
C
C   Now determine the N-X electron system information
C
C
         NBETA     = 29   ! Number of occupied beta  electrons (HOMO)
         NALPHA    = 29   ! Number of occupied alpha electrons (HOMO)
         BETA_ORB  = 1    ! MO number to swap to valence
         ALPHA_ORB = 0    ! MO number to swap to valence
         ASWAP  = .TRUE.  ! Swap an alpha MO
         BSWAP  = .TRUE.  ! Swap a  beta  MO
         BREAK_SYM = .TRUE.
         IF (ALPHA_ORB .EQ. 0) ASWAP = .FALSE.
         IF (BETA_ORB  .EQ. 0) BSWAP = .FALSE.
C
C   Symmetry breaking
C
         IF (BREAK_SYM) THEN

            WRITE (*,*) ' BREAKING SYMMETRY OF ORBITALS 1 AND 2 '

            DROOT2 = 1.D0 / DSQRT (2.0D0)

            DO MU = 1, ADIM1
               TMP (MU,1) = DROOT2 * (CA (MU,1) + CA (MU,2))
               TMP (MU,2) = DROOT2 * (CA (MU,1) - CA (MU,2))
            ENDDO

            DO MU = 1, ADIM1
               CA (MU,1) = TMP (MU,1)
               CA (MU,2) = TMP (MU,2)
            ENDDO

            DO MU = 1, ADIM1
               TMP (MU,1) = DROOT2 * (CB (MU,1) + CB (MU,2))
               TMP (MU,2) = DROOT2 * (CB (MU,1) - CB (MU,2))
            ENDDO

            DO MU = 1, ADIM1
               CB (MU,1) = TMP (MU,1)
               CB (MU,2) = TMP (MU,2)
            ENDDO

         ENDIF
C
C   Beta Swap
C
         IF (BSWAP) THEN

            WRITE (*,1000) ' SWAPPING BETA  ORBITAL ',BETA_ORB,
     +                     ' WITH ',NBETA
C
C
C   Put the orbital to be swapped in the temp virtual created spot
C
C
            DO MU = 1, BDIM1
               TMP (MU,NBETA) = CB (MU,BETA_ORB)
            ENDDO
C
C
C   Move every orbital above the ionized orbital down one
C
C
            IF (BETA_ORB .GT. 1) THEN
               DO Q = 1, BETA_ORB
                  DO MU = 1, BDIM1
                     TMP (MU,Q) = CB (MU,Q)
                  ENDDO
               ENDDO
            ENDIF

            DO Q = BETA_ORB, NBETA-1
               DO MU = 1, BDIM1
                  TMP (MU,Q) = CB (MU,Q+1)
               ENDDO
            ENDDO
C
C
C   Now copy the rearranged coefficient array in the the real CB
C
C
            DO Q = 1, NBETA
               DO MU = 1, BDIM1
                  CB (MU,Q) = TMP (MU,Q)
               ENDDO
            ENDDO

         ENDIF
C
C   Alpha Swap
C
         IF (ASWAP) THEN

            WRITE (*,1000) ' SWAPPING ALPHA ORBITAL ',ALPHA_ORB,
     +                     ' WITH ',NALPHA
C
C
C   Put the orbital to be swapped in the temp virtual created spot
C
C
            DO MU = 1, ADIM1
               TMP (MU,NALPHA) = CA (MU,ALPHA_ORB)
            ENDDO
C
C
C   Move every orbital above the ionized orbital down one
C
C
            IF (ALPHA_ORB .GT. 1) THEN
               DO P = 1, ALPHA_ORB
                  DO MU = 1, ADIM1
                     TMP (MU,P) = CA (MU,P)
                  ENDDO
               ENDDO
            ENDIF

            DO P = ALPHA_ORB, NALPHA-1
               DO MU = 1, ADIM1
                  TMP (MU,P) = CA (MU,P+1)
               ENDDO
            ENDDO
C
C
C   Now copy the rearranged coefficient array in the the real CB
C
C
            DO P = 1, NALPHA
               DO MU = 1, ADIM1
                  CA (MU,P) = TMP (MU,P)
               ENDDO
            ENDDO

         ENDIF
C
C
C   Break the symmetry in the approximate IP-VCCS version!
C
C
      ENDIF
      ELSE IF (read_type .EQ. 2) THEN

         BREAK_SYM = .TRUE.

         IF (BREAK_SYM) THEN

            IF (me .EQ. 0) 
     +         WRITE (*,*) ' BREAKING SYMMETRY OF ORBITALS 1 AND 2',
     +                     ' IN APPROXIMATE SCHEME '

            DROOT2 = 1.0D0 / DSQRT (2.0D0)

            DO MU = 1, ADIM1
               TMP (MU,1) = DROOT2 * (CA (MU,1) + CA (MU,2))
               TMP (MU,2) = DROOT2 * (CA (MU,1) - CA (MU,2))
            ENDDO

            DO MU = 1, ADIM1
               CA (MU,1) = TMP (MU,1)
               CA (MU,2) = TMP (MU,2)
            ENDDO

            DO MU = 1, BDIM1
               TMP (MU,1) = DROOT2 * (CB (MU,1) + CB (MU,2))
               TMP (MU,2) = DROOT2 * (CB (MU,1) - CB (MU,2))
            ENDDO

            DO MU = 1, BDIM1
               CB (MU,1) = TMP (MU,1)
               CB (MU,2) = TMP (MU,2)
            ENDDO

         ENDIF

      ENDIF ! (read_type)

 1000 FORMAT(A24,I3,A6,I3)
 1001 FORMAT (A4,I3,A6,F20.8)

      RETURN
      END
