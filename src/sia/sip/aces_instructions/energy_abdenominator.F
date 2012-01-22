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
      subroutine energy_abdenominator(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c   Divides each sample of the block given by the array argument by the 
c       eps = f_alpha(i)+f_beta(i)-f_alpha(a)-f_beta(a)
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'epsilon.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      if (nindex .eq. 2) then
         call epsab_divide2(x, bval(1),eval(1),type(1),
     *                    bdim(1), edim(1),
     *                    bval(2),eval(2),type(2),
     *                    bdim(2), edim(2),
     *                    epsilon,epsilonb)    
      else
         print *,'Error in energy_abdenominator: Result array ',
     *      'must have 2 indices.'
         print *,'array has ',nindex,' indices.'
         call abort_job()
      endif
      
      return
      end

      subroutine epsab_divide2(x, a1,a2,atype,na1,na2,i1,i2,itype,
     *                       ni1,ni2,epsilon, epsilonb)
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
            epsi = epsilon(i) + epsilonb(i)
      do a = a1,a2
            epsa = epsilon(a) + epsilonb(a)
         val =  x(a,i)
         eps = epsi - epsa
         if (eps .ne. 0.0) x(a,i) = val/eps

         if (dabs(eps) .lt. 1.0d-10) then 
            if (a .ne. i) write(6,*) ' Small denominator being
     *         eliminated in energy_denominator', 
     *         'a, i', a,i, 'den' ,eps, 'value', x(a,i) 
            x(a,i) = 0.0d0  
c           x(a,i) = val   
         endif 
c
      enddo
      enddo

      return
      end
