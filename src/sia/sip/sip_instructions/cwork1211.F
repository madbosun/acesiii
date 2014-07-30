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
      subroutine cwork1211(y,na1,na2,a1,a2,inda,
     *                     x1,nd1,nd2,ne1,ne2,
     *                     d1,d2,e1,e2, indx1,
     *                     x2,ni1,ni2,
     *                     i1,i2,indi, 
     *                     scr1, scr2,flopcount) 
c-------------------------------------------------------------------------
c   Performs a "1211" contraction: 
c      1 index output array  --> y 
c      2 index operand array --> x1 
c      1 index operand array --> x2 
c      1 index contraction.
c
c   I. e., contract the 2 indices of the 2nd operand array out of the 
c   first operand array into the 2-index output array.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      character*1 char(5)  
      character*2 char2a, char2b   
      character*5 char5 

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,i1,i2,j1,j2
      integer inda, indb, indc, indx1(2), indi, indj, h(3), perm(5) 
      integer flopcount

      double precision y(na1:na2)
      double precision x1(nd1:nd2,ne1:ne2)
      double precision x2(ni1:ni2)
      double precision scr1(*), scr2(*), x  

      integer ia, ib, ii, ij, ix(2)
      integer i, j, a, b, c, d, e, f 
#ifdef __crayx1
      integer*8 ncount, m, n, k, lda, ldb, ldc 
#else
      integer ncount, m, n, k, lda, ldb, ldc 
#endif

      logical check
      check = .false.

c      write(6,*) ' Entering cwork1211 ', inda, indx1(1), indx1(2)  

        do d=nd1,nd2
        y(d) = 0.d0
        enddo 


      if (indx1(1) .eq. inda) then 
         do e = ne1, ne2  
         do d = nd1, nd2  
            y(d) = y(d) + x1(d,e)*x2(e) 
         enddo 
         enddo 
      endif 

      if (indx1(2) .eq. inda) then 
         do d = nd1, nd2  
         do e = ne1, ne2  
            y(d) = y(d) + x1(e,d)*x2(e) 
         enddo 
         enddo 
      endif 

      return
      end

