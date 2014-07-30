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
      subroutine twork2111(y,na1,na2,nb1,nb2,
     *                     a1,a2,b1,b2,indy,
     *                     x1,ne1,ne2,
     *                     e1,e2,indx1,
     *                     x2,ni1,ni2,
     *                     i1,i2,indx2,
     *                     flopcount) 
c-------------------------------------------------------------------------
c   Performs a "4222" tensor contraction: 
c      4 index output array 
c      2 index operand array
c      2 index operand array
c      2 index contraction.
c
c--------------------------------------------------------------------------
      implicit none

      integer na1,na2,nb1,nb2,ne1,ne2,
     *        ni1,ni2 
      integer a1,a2,b1,b2,e1,e2,i1,i2
      integer indy(2), indx1(1), indx2(1)
      integer flopcount

      double precision y(na1:na2,nb1:nb2)
      double precision x1(ne1:ne2)
      double precision x2(ni1:ni2)

      integer i, j, a, b, c, d

      flopcount = (a2-a1+1)*(b2-b1+1)
     
c----------------------------------------------------------------------------
c   ia is the 1st index of x1.
c----------------------------------------------------------------------------

      if (indx1(1) .eq. indy(1)) then 
         do b = b1,b2
         do a = a1,a2
            y(a,b) = x1(a)*x2(b)
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. indy(2)) then 
         do b = b1,b2
         do a = a1,a2
            y(a,b) = x1(b)*x2(a)
         enddo 
         enddo 

         return
      endif

      print *,'twork2111: Your tensor operation does not fit',
     *    ' any valid pattern'
      call abort_job() 

      return
      end

