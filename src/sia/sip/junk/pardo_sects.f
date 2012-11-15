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
      subroutine pardo_sects(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute pardo_sects scalar1 scalar2
c          
c   scalar1 is type of calculation.
c   scalar2 is output.
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'

      double precision x
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      integer i, j, k
      integer vcount
      integer oseg, vseg, nsects

      if (nindex .ne. 0) then
         print *,'Error: pardo_sects at line ',current_line
         print *,'First arg. must be a scalar.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: pardo_sects at line ',current_line,' not ',
     *           'called with scalar in 2nd arg.'
         call abort_job()      
      endif

c --------------------------------------------------------------------------- 
c Define the output scalar --> number of segmentation 
c --------------------------------------------------------------------------- 

      oseg = x2  ! The number of occupied triplets 
      vseg = nalpha_virtual/sip_mx_virt_segsize
c
c --------------------------------------------------------------------------- 
c RHF AAA triples calculation 
c --------------------------------------------------------------------------- 
c 
      if (x .eq. 1) then

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = i, vseg ! nalpha_virtual 
         do k = j, vseg ! nalpha_virtual 
            vcount = vcount + 1
         enddo
         enddo
         enddo

         nsects = my_company_size/vcount + 1
         if (oseg .lt. nsects) nsects = oseg
         if (nsects .gt. 14.0) nsects = 14.0
         if (nsects .lt. 4.0) nsects = 4.0
         x2 = nsects  

         return

      endif
c
c --------------------------------------------------------------------------- 
c RHF AAB triples calculation 
c --------------------------------------------------------------------------- 
c 
      if (x .eq. 2) then

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = i, vseg ! nalpha_virtual 
         do k = 1, vseg ! nalpha_virtual 
            vcount = vcount + 1
         enddo 
         enddo 
         enddo 

         nsects = my_company_size/vcount + 1
         if (oseg .lt. nsects) nsects = oseg
         if (nsects .gt. 14.0) nsects = 14.0
         if (nsects .lt. 4.0) nsects = 4.0
         x2 = nsects  

         return 
      
      endif
c
c --------------------------------------------------------------------------- 
c UHF aaa ring calculation 
c --------------------------------------------------------------------------- 
c 
      if (x .eq. 3) then

         vseg = nalpha_virtual/sip_mx_virt_segsize
         oseg = nalpha_occupied/sip_mx_occ_segsize

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = 1, oseg ! nalpha_occupied  
            vcount = vcount + 1
         enddo 
         enddo 

         nsects = my_company_size/vcount + 1
         if (nsects .gt. 5.0) nsects = 5.0
         x2 = nsects  

         return 
      
      endif
c
c ---------------------------------------------------------------------------
c UHF aaa ring calculation
c ---------------------------------------------------------------------------
c
      if (x .eq. 3) then

         vseg = nalpha_virtual/sip_mx_virt_segsize
         oseg = nalpha_occupied/sip_mx_occ_segsize

         if (vseg .lt. 1) vseg = 1
         if (oseg .lt. 1) oseg = 1

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = 1, oseg ! nalpha_occupied
            vcount = vcount + 1
         enddo
         enddo

         nsects = my_company_size/vcount + 1
         if (nsects .gt. 5.0) nsects = 5.0
         x2 = nsects

         return

      endif

      
      
      call c_flush_stdout()
      return
      end        
         

