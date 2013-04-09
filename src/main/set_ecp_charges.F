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

      subroutine set_ecp_charges(max_centers,natoms,icrcor)
      implicit none

#include "mxatms.par"

      integer natoms  
      integer max_centers 
      integer charge(max_centers) 
      integer icrcor(mxatms) 
      integer nlines, ios, i  
      character*80 line
      character*2 name(max_centers) 

c Open the ZMAT file 

      open (unit=23, file='ZMAT', status='OLD',
     *      err = 200, iostat = ios)

  200 continue 
      if (ios .ne. 0) then
         print *,'Error: Cannot open ZMAT file.'
         call abort_job()
      endif


c Read atom data from the ZMAT 

      natoms = 0 
      nlines = 0 
  300 continue
      line = ' '
      read (23,9000, end = 400) line
      nlines = nlines + 1 
      if (nlines .gt. 1) then 
         if (line .eq. '  ') go to 400 
         natoms = natoms + 1 
111   continue 
      endif 

      go to 300

  400 continue 

      close(23) 

#ifdef _DEBUG_LVL0
      write(6,*) '@set_ecp_charges: Number of ECP ATOMS = ', natoms
#endif

9000  FORMAT(a80) 
      return
      end

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------

