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
      subroutine twork642(y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                     ne1,ne2,nf1,nf2, 
     *                     a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,
     *                     f1,f2,indy,
     * 
     *                     x1,nk1,nk2,nl1,nl2,nm1,nm2,nn1,nn2,
     *                     k1,k2,l1,l2,m1,m2,n1,n2, indx1,
     * 
     *                     x2,ni1,ni2,nj1,nj2,
     *                     i1,i2,j1,j2,indx2,
     *                     flopcount) 
c-------------------------------------------------------------------------
c   Performs a "642" tensor contraction: 
c      6 index output array 
c      4 index operand array
c      2 index operand array
c
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'parallel_info.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        nk1,nk2,nl1,nl2,nm1,nm2,nn1,nn2,ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,k1,k2,l1,
     *        l2,m1,m2,n1,n2,i1,i2,j1,j2
      integer indy(6), indx1(4), indx2(2)
      integer flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,
     *                   ne1:ne2,nf1:nf2) 
      double precision x1(nk1:nk2,nl1:nl2,nm1:nm2,nn1:nn2)
      double precision x2(ni1:ni2,nj1:nj2)

      integer ia, ib, ic, id, ie, if, ia2, ib2, ic2, id2, 
     *        ie1, ie2, if1, if2
      integer ix(0:2), ix2(0:2)
      integer i, j, k, l, m, n, a, b, c, d, e, f 

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y
c   indices.
c---------------------------------------------------------------------------

      do i = 0,2 
         ix(i) = 0
         ix2(i) = 0
      enddo

      ia = 0
      ib = 0
      ic = 0
      id = 0
      ie = 0
      if = 0
      do j = 1, 4
      do i = 1, 6
         if (indx1(j) .eq. indy(i)) then
            if (i .eq. 1) then
               ia = indy(1)
            else if (i .eq. 2) then
               ib = indy(2)
            else if (i .eq. 3) then
               ic = indy(3)
            else if (i .eq. 4) then
               id = indy(4)
            else if (i .eq. 5) then
               ie = indy(5)
            else 
               if = indy(6)
            endif
         endif 
      enddo
      enddo

c---------------------------------------------------------------------------
c   Find which indices of the "x2" operand match the various y
c   indices.
c---------------------------------------------------------------------------

      ia2 = 0
      ib2 = 0
      ic2 = 0
      id2 = 0
      ie2 = 0
      if2 = 0
      do j = 1, 2
      do i = 1, 6
         if (indx2(j) .eq. indy(i)) then
            if (i .eq. 1) then
               ia2 = indy(1)
            else if (i .eq. 2) then
               ib2 = indy(2)
            else if (i .eq. 3) then
               ic2 = indy(3)
            else if (i .eq. 4) then
               id2 = indy(4)
            else if (i .eq. 5) then
               ie2 = indy(5)
            else 
               if2 = indy(6)
            endif
         endif 
      enddo
      enddo

      flopcount = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)*(d2-d1+1)
     *          * (e2-e1+1)*(f2-f1+1) 

c---------------------------------------------------------------------------
c Permutation (abcdef) = (abcd)*(ef) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(1) .and. indx1(2) .eq. indy(2) .and.
     *    indx1(3) .eq. indy(3) .and. indx1(4) .eq. indy(4) .and. 
     *    indx2(1) .eq. indy(5) .and. indx2(2) .eq. indy(6)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,b,c,d)*x2(e,f)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (abcd)*(fe) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(1)  .and. indx1(2) .eq. indy(2) .and.
     *    indx1(3) .eq. indy(3)  .and. indx1(4) .eq. indy(4) .and. 
     *    indx2(1) .eq. indy(6) .and. indx2(2) .eq. indy(5)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,b,c,d)*x2(f,e)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cbef)*(ad) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(2) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(1) .and. indx2(2) .eq. indy(4)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,b,e,f)*x2(a,d)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (adef)*(cb) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(1)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(3) .and. indx2(2) .eq. indy(2)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,d,e,f)*x2(c,b)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cdef)*(ab) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(1) .and. indx2(2) .eq. indy(2)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,d,e,f)*x2(a,b)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cbef)*(ad) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(2) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(1) .and. indx2(2) .eq. indy(4)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,b,e,f)*x2(a,d)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cdef)*(ab) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(1) .and. indx2(2) .eq. indy(2)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,d,e,f)*x2(a,b)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cdeb)*(af) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(2) .and. 
     *    indx2(1) .eq. indy(1) .and. indx2(2) .eq. indy(6)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,d,e,b)*x2(a,f)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cdaf)*(eb) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(1)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(5) .and. indx2(2) .eq. indy(2)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,d,a,f)*x2(e,b)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (abef)*(cd) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(1)  .and. indx1(2) .eq. indy(2) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(3) .and. indx2(2) .eq. indy(4)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,b,e,f)*x2(c,d)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cbaf)*(ed) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(2) .and.
     *    indx1(3) .eq. indy(1)  .and. indx1(4) .eq. indy(6) .and. 
     *    indx2(1) .eq. indy(5) .and. indx2(2) .eq. indy(4)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,b,a,f)*x2(e,d)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (adeb)*(cf) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(1)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(5)  .and. indx1(4) .eq. indy(2) .and. 
     *    indx2(1) .eq. indy(3) .and. indx2(2) .eq. indy(6)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,d,e,b)*x2(c,f)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

c---------------------------------------------------------------------------
c Permutation (abcdef) = (cdab)*(ef) 
c---------------------------------------------------------------------------
     
      if (indx1(1) .eq. indy(3)  .and. indx1(2) .eq. indy(4) .and.
     *    indx1(3) .eq. indy(1)  .and. indx1(4) .eq. indy(2) .and. 
     *    indx2(1) .eq. indy(5) .and. indx2(2) .eq. indy(6)) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(c,d,a,b)*x2(e,f)
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

C
C totally unoptimized crap here
C
         if (me.eq.0) then
         write(6,*) 'Using unoptimized cwork642 at line number ',
     *            current_line
         endif
         do f = nj1,nj2
         do e = ni1,ni2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(b,c,d,e)*x2(a,f)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         return


c---------------------------------------------------------------------------

      print *,'twork642: Your tensor operation does not fit',
     *    ' any valid pattern'
      write(6,*) ' Line number :', current_line 
      call abort_job() 

      return
      end

