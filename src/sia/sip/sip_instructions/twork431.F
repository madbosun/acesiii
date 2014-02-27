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
c     subroutine twork431(y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
c    *                     a1,a2,b1,b2,c1,c2,d1,d2,indy,
c    *                     x1,ne1,ne2,nf1,nf2,ng1,ng2,
c    *                     e1,e2,f1,f2,g1,g2,indx1,
c    *                     x2,ni1,ni2,
c    *                     i1,i2,indx2,
c    *                     flopcount) 
      subroutine twork431(y,ne1,ne2,nf1,nf2,ng1,ng2,ni1,ni2,
     *                     e1,e2,f1,f2,g1,g2,i1,i2,indy,
     * 
     *                     x1,na1,na2,nb1,nb2,nc1,nc2,
     *                     a1,a2,b1,b2,c1,c2,indx1,
     * 
     *                     x2,nd1,nd2,
     *                     d1,d2,indx2,
     *                     flopcount) 
c-------------------------------------------------------------------------
c   Performs a "431" tensor contraction: 
c      4 index output array 
c      3 index operand array
c      1 index operand array
c
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,ng1,ng2,
     *        ni1,ni2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,i1,i2
      integer indy(4), indx1(3), indx2(1)
      integer flopcount

      double precision y(ne1:ne2,nf1:nf2,ng1:ng2,ni1:ni2)
      double precision x1(na1:na2,nb1:nb2,nc1:nc2)
      double precision x2(nd1:nd2)

      integer ia, ib, ic, id, ia2, ib2, ic2, id2
      integer i, j, a, b, c, d
      integer perm(4), h(4), ncount   

      character*4 char4 
      character*1 char(4) 

      flopcount = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)*(d2-d1+1)
     
c---------------------------------------------------------------------------
c   Determine permutation required to put output array into desired order 
c---------------------------------------------------------------------------

      char(1) = '1'
      char(2) = '2'
      char(3) = '3'
      char(4) = '4'

      h(1) = indx1(1) 
      h(2) = indx1(2)  
      h(3) = indx1(3)  
      h(4) = indx2(1)  

      ncount = 0
      do i = 1, 4
      do j = 1, 4
         if (h(i) .eq. indy(j)) then 
            perm(j) = i
            ncount = ncount + 1
         endif
      enddo
      enddo

      if (ncount .ne. 4) then
         write(6,*) ' Problem with twork431 at line :', current_line
         call abort_job()
      endif 


      char4 = char(perm(1))//char(perm(2))//char(perm(3))//char(perm(4))

c
      if (char4 .eq. '1234') then
c          write(6,*) '11111111111' 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
      if (char4 .eq. '1243') then
c           write(6,*) '22222222222' 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,d,c) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
      if (char4 .eq. '2143') then 
c           write(6,*) '33333333333'
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(b,a,d,c) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c

      if (char4 .eq. '1324') then 
c           write(6,*) '4444444444'
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,c,b,d) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
      if (char4 .eq. '3214') then
c           write(6,*) '5555555555' 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(c,b,a,d) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
c-----------------------

      if (char4 .eq. '4231') then
c        write(6,*) '++++++++++++++++++++++++++++'  

c        write(6,*) ' ne nd:', ne1, ne2, ':', nd1, nd2 
c        write(6,*) ' e d:', e1, e2, ':', d1, d2
c        write(6,*) ' Problem with twork431 at line :', current_line
       
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(d,b,c,a) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c

      if (char4 .eq. '2134') then 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(b,a,c,d) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c

c permutation 3142 
      if (char4 .eq. '3142') then 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(c,a,d,b) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
c permutation 4231  
      if (char4 .eq. '4231') then 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(d,b,c,a) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
c permutation 3124 
      if (char4 .eq. '3124') then 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(c,a,b,d) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
c permutation 2341 
      if (char4 .eq. '2341') then 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(b,c,d,a) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c
c permutation 4123  
      if (char4 .eq. '4123') then 
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(d,a,b,c) = x1(a,b,c)*x2(d)
         enddo 
         enddo 
         enddo 
         enddo 
         return
      endif 
c

      print *,'twork431: Your tensor operation does not fit',
     *    ' any valid pattern; This sucks.', char4 
      call abort_job() 

      return
      end

