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
      subroutine c1_full_print(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Prints C1 a and C1b vectors
c   execute c1a c1b
c 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'

      double precision x(*) 
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      integer i, n, i1, i2, a1, a2, j1, j2, b1, b2
 
      i1 = bdim(1)
      i2 = edim(1)
      a1 = bdim(2)
      a2 = edim(2)
      j1 = bdim2(1)
      j2 = edim2(1)
      b1 = bdim2(2)
      b2 = edim2(2)

      call uuu (x, x2, i1, i2, a1, a2, j1, j2, b1, b2)

      return
      end




      subroutine uuu (c1a, c1b, i1, i2, a1, a2, j1, j2, b1, b2)

      implicit none
      include 'mpif.h'
      include 'parallel_info.h'

      double precision stemp, sumb
      integer i, n, a, b, i1, i2, a1, a2, j1, j2, b1, b2
      double precision c1a (i1:i2,a1:a2), c1b (j1:j2,b1:b2)

10    FORMAT(2x,'Singles amplitudes above 0.1')
20    FORMAT(2x,60('-'))
30    FORMAT(2x,I8,x,I8,x,F16.8)
40    FORMAT(2x,'A',10x,'I',10x,'C1-Alpha (A,I)')
50    FORMAT(2x,'A',10x,'I',10x,'C1-Beta  (A,I)')

       open (unit=66,file="summary.out")
       if (me .eq. my_company_rank)then

         write(66,*) ' '
         write(66,10)
         write(66,20)
         write(66,40)
c
         do b = a1, a2
         do a = i1, i2
         stemp=dabs(c1a(a,b))
c
         if(stemp .gt. 1.0d-1)then
           write(66,30) a, b, c1a(a,b)
         endif

         enddo
         enddo

         write(66,20)
         sumb = 0.0d0
         do b = b1, b2
         do a = j1, j2
            sumb= sumb+dabs(c1b(a,b))
         enddo
         enddo 

         if (sumb .gt. 0.0d0) then
            write(66,50)

            do b = b1, b2
            do a = j1, j2
               stemp=dabs(c1b(a,b))
c
               if(stemp .gt. 1.0d-1)then
                  write(66,30) a, b, c1b(a,b)
               endif

            enddo
            enddo
         endif
         write(66,20)

c
       endif

       return
       end
