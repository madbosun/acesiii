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
      subroutine cwork2311(y,na1,na2,nb1,nb2,a1,a2,b1,b2,inda, indb, 
     *                     x1,nc1,nc2,nd1,nd2,ne1,ne2,
     *                     c1,c2,d1,d2,e1,e2, indx1,
     *                     x2,ni1,ni2,
     *                     i1,i2,indi, 
     *                     scr1, scr2,flopcount) 
c-------------------------------------------------------------------------
c   Performs a "2311" contraction: 
c      2 index output array
c      3 index operand array
c      1 index operand array
c      1 index contraction.
c
c   I. e., contract the 2 indices of the 2nd operand array out of the 
c   first operand array into the 2-index output array.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,i1,i2,j1,j2
      integer inda, indb, indx1(3), indi, indj
      integer flopcount

      double precision y(na1:na2,nb1:nb2)
      double precision x1(nc1:nc2,nd1:nd2,ne1:ne2)
      double precision x2(ni1:ni2)
      double precision scr1(*), scr2(*) 

      character*2 char2
      character*1 char(3) 
      integer h(3), perm(3), indy(2)  
      double precision etemp 

      integer ia, ib, ii, ij, ix(4)
      integer i, j, a, b, c, d, e, f 
#ifdef __crayx1
      integer*8 ncount, m, n, k, lda, ldb, ldc 
#else
      integer ncount, m, n, k, lda, ldb, ldc 
#endif

      logical check
      check = .false.

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y and x2 
c   indices.
c---------------------------------------------------------------------------

c      flopcount = 2 * (i2-i1+1)*(j2-j1+1)*(b2-b1+1)*(a2-a1+1)
 
       flopcount = 2 * (i2-i1+1)*(b2-b1+1)*(a2-a1+1)

c---------------------------------------------------------------------------
c   Determine permutation required to put output array into desired order 
c---------------------------------------------------------------------------

      indy(1) = inda 
      indy(2) = indb 

      char(1) = '1'
      char(2) = '2'
      char(3) = '3'

      h(1) = indx1(1)
      h(2) = indx1(2)
      h(3) = indx1(3)

      do i = 1, 3
      do j = 1, 2
         if (h(i) .eq. indy(j)) perm(j) = i
      enddo
      enddo

      char2 = char(perm(1))//char(perm(2))

c---------------------------------------------------------------------------
c   Only do ONE permutation for now  
c---------------------------------------------------------------------------

c      if (char2 .eq. '12') then 
         DO a = a1, a2       
         DO b = b1, b2       
            etemp = 0.d0 
            DO c = c1, c2 
               etemp = etemp + x1(a,b,c)*x2(c)  
            ENDDO 
            y(a,b) = etemp 
         ENDDO 
         ENDDO 
         return 
c      endif  


c       write(6,*) ' PERMUTATION ', char4, 
c     *     ' NOT OPTIMIZED IN cwork2422 line ',current_line 

      return
      end

