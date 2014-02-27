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
      subroutine cwork3411(y,na1,na2,nb1,nb2,nc1,nc2,
     *                     a1,a2,b1,b2,c1,c2,inda, indb,
     *                     indc, 
     * 
     *                     x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *                     e1,e2,f1,f2,g1,g2,h1,h2,
     *                     indx1, 
     * 
     *                     x2,nd1,nd2,
     *                     d1,d2,indd, 
     *                     cind, flopcount,scr1, scr2)
c-------------------------------------------------------------------------
c   Performs a "3411" contraction:
c      3 index output array
c      4 index operand array
c      1 index operand array
c      1 index contraction.
c
c   I. e., contract one of the the 2 indices of the 2nd operand array out
c   of the first operand array, replacing the index by the 2nd operand's
c   2nd index.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2 
      integer ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2
      integer ni1, ni2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,
     *        i1,i2
      integer inda, indb, indc, indx1(4), indd 
      integer cind(2), flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2)
      double precision x2(nd1:nd2)

      double precision scr1(*)
      double precision scr2(*)
      integer next

      integer ia, ib, ic, icontract, ix(3), iy
      integer i, j, a, b, c, d 
      integer indy(3), map(3)
      integer ie, if, ig
      integer e, f, g 
      integer alast, dlast
      integer iprint
      double precision cval

      character*1 char(4) 
      character*4 char4 
      integer h(4), perm(4) 

      integer*8 ixx, c_loc64

#ifdef __crayx1
      integer*8 m, n, k, lda, ldb, ldc
#else
      integer m, n, k, lda, ldb, ldc
#endif

      char(1) = '1'
      char(2) = '2'
      char(3) = '3'
      char(4) = '4'

      h(1) = inda
      h(2) = indb
      h(3) = indc
      h(4) = indd

c---------------------------------------------------------------------------
c   Determine permutation of array 1 
c---------------------------------------------------------------------------



c---------Initialization of array-------------------------------------------
      do a = a1, a2
      do b = b1, b2
      do c = c1, c2
      Y(a,b,c) = 0.d0
      enddo
      enddo
      enddo

c--------------------------------------------------------------------------
      do i = 1, 4
      do j = 1, 4
         if (h(i) .eq. indx1(j)) perm(j) = i
      enddo
      enddo

      char4 = char(perm(1))//char(perm(2))//char(perm(3))//char(perm(4)) 

c---------------------------------------------------------------------------
c   Y(a,b,c) = X1(a,b,c,d)*x2(d) 
c---------------------------------------------------------------------------

      if (char4 .eq. '1234') then 

         do d = d1, d2 
         do c = c1, c2 
         do b = b1, b2 
         do a = a1, a2 
            Y(a,b,c) = Y(a,b,c) + x1(a,b,c,d)*x2(d) 
         enddo 
         enddo 
         enddo 
         enddo 
         return 

      endif 

c---------------------------------------------------------------------------
c   Y(a,b,c) = X1(a,b,d,c)*x2(d) 
c---------------------------------------------------------------------------

      if (char4 .eq. '1243') then 

         do d = d1, d2 
         do c = c1, c2 
         do b = b1, b2 
         do a = a1, a2 
            Y(a,b,c) = Y(a,b,c) + x1(a,b,d,c)*x2(d) 
         enddo 
         enddo 
         enddo 
         enddo 
         return 

      endif 

c---------------------------------------------------------------------------
c   Y(a,b,c) = X1(a,d,b,c)*x2(d) 
c---------------------------------------------------------------------------

      if (char4 .eq. '1423') then 

         do d = d1, d2 
         do c = c1, c2 
         do b = b1, b2 
         do a = a1, a2 
            Y(a,b,c) = Y(a,b,c) + x1(a,d,b,c)*x2(d) 
         enddo 
         enddo 
         enddo 
         enddo 
         return 

      endif 

c---------------------------------------------------------------------------
c   Y(a,b,c) = X1(d,a,b,c)*x2(d) 
c---------------------------------------------------------------------------

      if (char4 .eq. '4123') then 

         do c = c1, c2 
         do b = b1, b2 
         do a = a1, a2 
         do d = d1, d2 
            Y(a,b,c) = Y(a,b,c) + x1(d,a,b,c)*x2(d) 
         enddo 
         enddo 
         enddo 
         enddo 
         return 

      endif 




      print *,'Error: Invalid contraction at line ',current_line
      write(6,*) 'char4 = ', char4
      call abort_job()  
      return
      end
       
