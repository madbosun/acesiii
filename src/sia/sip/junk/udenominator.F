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
      subroutine udenominator(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c   Divides each sample of the block given by the array argument by the
c   denominator
c       eps = epsilon(p)-epsilon(q)
c   The p=q value is set to zero.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      if (nindex .eq. 4) then
         call epsu_divide(x, bval(1),eval(1),type(1),
     *                    bval(2),eval(2),type(2),
     *                    bval(3),eval(3),type(3),
     *                    bval(4),eval(4),type(4), 
     *                    epsilon,epsilonb)
      else if (nindex .eq. 2) then
         call epsu_divide2(x, bval(1),eval(1),type(1),
     *                     bdim(1), edim(1),
     *                     bval(2),eval(2),type(2),
     *                     bdim(2), edim(2),
     *                     epsilon,epsilonb)    
      else
         print *,'Error in udenominator: Result array ',
     *      'must have either 4 or 2 indices.'
         print *,'array has ',nindex,' indices.'
         call abort_job()
      endif
      
      return
      end

      subroutine epsu_divide(x, a1,a2,atype,
     *                       i1,i2,itype,b1,b2,btype,
     *                       j1,j2, jtype, epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j
      double precision eps, epsa, epsb, epsi, epsj
      double precision val

      do j = j1,j2
         if (jtype .eq. mobindex) then
            epsj = epsilonb(j)
         else
            epsj = epsilon(j)
         endif
      do b = b1,b2
         if (btype .eq. mobindex) then
            epsb = epsilonb(b)
         else
            epsb = epsilon(b)
         endif
      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)
         else
            epsi = epsilon(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif
         val =  x(a,i,b,j)
         eps = epsi + epsj - epsa - epsb
         if (eps .ne. 0.0) x(a,i,b,j) = val/eps

c A more robust check may be needed here VFL 

      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine epsu_divide2(x, a1,a2,atype,na1,na2,i1,i2,itype,
     *                        ni1,ni2,epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,i1,i2
      integer atype, itype
      integer na1, na2, ni1, ni2
      double precision x(na1:na2,ni1:ni2)

      double precision epsilon(*), epsilonb(*)

      integer a,i
      double precision eps, epsa, epsi
      double precision val

      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)
         else
            epsi = epsilon(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif
         val =  x(a,i)
         eps = epsi - epsa

c        if (eps .ne. 0.0) x(a,i) = val/eps

         if (dabs(eps) .gt. 1.0d-10) x(a,i) = val/eps  
         if (dabs(eps) .lt. 1.0d-10) then 
            x(a,i) = 0.0d0  
            if (a .ne. i) write(6,*) ' Small denominator may need  
     *         to be eliminated in energy_denominator', 
     *         'a, i', a,i, 'den' ,eps, 'value', x(a,i) 
            x(a,i) = 0.0d0  
         endif 

c If indices are the same do not modify array 

         if (a .eq. i) x(a,i) = 0.0d0  
c
      enddo
      enddo

      return
      end
