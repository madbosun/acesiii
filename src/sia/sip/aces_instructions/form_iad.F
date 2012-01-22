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
      subroutine form_iad(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      integer na1, na2, nb1, nb2

      na1 = edim(1)-bdim(1) + 1
      na2 = edim(2)-bdim(2) + 1
      nb1 = edim2(1)-bdim2(1) + 1
      nb2 = edim2(2)-bdim2(2) + 1
      call do_formiad(x,x2,nindex,na1,na2,nb1,nb2)

      return 
      end 

      subroutine do_formiad(array1,array2,nindex,na1,na2,nb1,nb2) 
      implicit none

      include 'int_gen_parms.h'

c----------------------------------------------------------------------------
c
c----------------------------------------------------------------------------

      integer nindex, na1, na2, nb1, nb2  
      double precision array1(na1,na2)
      double precision array2(nb1,nb2)
      double precision etemp  
      integer a, b, p, n1, n2, noff   

      n1 = min0(na1,nb1) 
      n2 = min0(na2,nb2) 

c     write(6,*) ' RANGE OF INDICES IN FORM_AID ' 
c     write(6,*) ' Array1 :', na1, na2   
c     write(6,*) ' Array2 :', nb1, nb2   

      do a = 1, n2   
      do b = 1, n2   
            etemp = 0.0 
            do p = 1, nalpha_occupied + nbeta_occupied  
               etemp = etemp + array1(a,p)*array1(b,p) 
            enddo
            array2(a,b) = 1.0*etemp  
c     write(6,*) ' a, b', a, b, array2(a,b) 
      enddo
      enddo

      return
      end

