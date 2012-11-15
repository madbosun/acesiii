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
      subroutine nuc_dipole_moment(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: nuc_dipole_moment scalar1 scalar2
c   scalar1 indicates x, y, or z moment to be computed.
c   scalar2 is the result.
c--------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'
      include 'interpreter.h'
      include 'trace.h'

	double precision x 
	integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)  

      double precision DNUC,XCOORD,YCOORD,ZCOORD

      integer imoment, natoms, i

      if (dryrun) return

      if (nindex .ne. 0) then
         print *,'Error: nuc_dipole_moment at line ',current_line
         print *,'First arg. must be a scalar.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: nuc_dipole_moment at line ',current_line,
     *           'must be called with scalar in 2nd arg.'
         call abort_job()
      endif

C
C
C   Compute the nuclear contribution to the dipole moment
C
C 
      if ( abs (x - 1) .LT. 0.0001 ) imoment = 1
      if ( abs (x - 2) .LT. 0.0001 ) imoment = 2
      if ( abs (x - 3) .LT. 0.0001 ) imoment = 3

      natoms = x2
      x2 = 0.0D0
      DNUC = 0.0D0
      IF (imoment .EQ. 1) THEN
         DO I = 1, natoms
            XCOORD = acenter(I,1)
            DNUC = DNUC + charge(I)*XCOORD
         ENDDO ! I
      ELSE IF (imoment .EQ. 2) THEN
         DO I = 1, natoms
            YCOORD = acenter(I,2)
            DNUC = DNUC + charge(I)*YCOORD
         ENDDO ! I
      ELSE IF (imoment .EQ. 3) THEN
         DO I = 1, natoms
            ZCOORD = acenter(I,3)
            DNUC = DNUC + charge(I)*ZCOORD
         ENDDO ! I
      ELSE
         WRITE(*,*) "Wrong type of moment; should be 1, 2, or 3."
         call abort_job()
      ENDIF
      x2 = DNUC

      return 
      end 


