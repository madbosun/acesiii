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
      subroutine copy_ff(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c
c   The elements of Array1 are put into Array2 for overlaping indeces. 
c
c   The instruction is executed as follows:
c   
c   execute copy_ff Array1 Array2
c
c   where Array1 is the full input matrix. 
c         Array2 is the same ar Array1 where indeces overlap.  
c
c   Both Array1 and Array2 must be "static" arrays.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      if (nindex .eq. 2) call copy2ff(x,x2,nindex,
     *                              bdim(1),edim(2),bdim(1),edim(2))
      return
      end

      subroutine copy2ff(array1,array2,nindex,na1,na2,nb1,nb2) 
      implicit none

      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'trace.h'

c----------------------------------------------------------------------------
c
c The array1 is assumed to be of type alpha-alpha (VaD) or (SD). The elements 
c of the array1 are put into array2 with the FIRST index offset by the
c number of singly occupied orbitals. Array2 is of type beta-beta (VbD)  
c 
c----------------------------------------------------------------------------

      integer nindex, na1, na2, nb1, nb2  
      double precision array1(na1,na2)
      double precision array2(nb1,nb2)
      integer a, b, n1, n2, noff   

      n1 = max0(na1,nb1) 
      n2 = max0(na2,nb2) 

c      write(6,*) ' RANGE OF INDICES IN COPY2 PROC ',me,' LINE ',
c     *    current_line 
c      write(6,*) ' Array1 :', na1, na2   
c      write(6,*) ' Array2 :', nb1, nb2   

      noff = 0 ! This routine stricktly copies one array into another 

      if (na1 .ne. nb1) then 
          write(6,*) ' Error in copy_ff. Arrays have different sizes ',
     *                 na1, nb1  
         call abort_job()
      endif 

      if (na2 .ne. nb2) then 
          write(6,*) ' Error in copy_ff. Arrays have different sizes ',
     *                 na2, nb2 
         call abort_job()
      endif 

      if (n1 .ne. n2) then 
          write(6,*) ' Error in copy_ff. Arrays have different sizes ',
     *                 n1, n2 
         call abort_job()
      endif 

      do b = 1, n2  
         do a = 1, n1  
            array2(a,b) = array1(a,b)
c           write(6,*) ' COPYFF:', a+noff, b, array1(a,b), 
c    *                                        array2(a,b) 
         enddo
      enddo

100   continue 

      return
      end

