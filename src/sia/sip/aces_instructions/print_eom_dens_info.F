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
      subroutine print_eom_dens_info(x,nindex, type, bval,
     *                              eval, bdim, edim,
     *                      x2, nindex2, type2, bval2,
     *                              eval2, bdim2, edim2)
c--------------------------------------------------------------------------
c   Usage: execute print_eom_dens_info array scalar
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

      if (nindex .ne. 2) then
         print *,'Error: print_eom_dens_info at line ',current_line
         print *,'First arg. must be a block of a 2-index array.'
         call abort_job()
      endif

      if (nindex2 .ne. 0) then
         print *,'Error: print_eom_dens_info at line ',current_line,
     *           'must be called with scalar in 2nd arg.'
         call abort_job()      
      endif

      call doprint_eom_dens_info(x, x2,
     *             bval(1),eval(1),bval(2),eval(2))

      return 
      end 

      subroutine doprint_eom_dens_info(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   

      if (me .eq. 0)then
         if (sval .eq. 1) then 
            write(6,*) ' ' 
            write(6,99)
            write(6,*) ' '
            write(6,98) 'Relevant data for this root '  
            write(6,99)
            write(6,*) ' '
            do a = a1, a2
               stemp=array1(a,a)
               write(6,100) stemp, stemp*27.2113957
            enddo
         endif 

         if (sval .eq. 2) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,97) stemp
            enddo
         endif

         if (sval .eq. 3) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                              X-DIRECTION '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Right",stemp
            enddo
         endif

         if (sval .eq. 4) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Left ",stemp
            enddo
         endif

         if (sval .eq. 5) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,102) stemp
            enddo
         endif

         if (sval .eq. 6) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,103) stemp
            enddo
         endif

         if (sval .eq. 7) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                              Y-DIRECTION '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Right",stemp
            enddo
         endif

         if (sval .eq. 8) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Left ",stemp
            enddo
         endif

         if (sval .eq. 9) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,102) stemp
            enddo
         endif

         if (sval .eq. 10) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,103) stemp
            enddo
         endif

         if (sval .eq. 11) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                              Z-DIRECTION '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Right",stemp
            enddo
         endif

         if (sval .eq. 12) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Left ",stemp
            enddo
         endif

         if (sval .eq. 13) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,102) stemp
            enddo
         endif

         if (sval .eq. 14) then 
            do a = a1, a2
               stemp=array1(a,a)
               write(6,103) stemp
            enddo
         endif 

         if (sval .eq. 15) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               XX-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 16) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 17) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               YY-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 18) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 19) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               ZZ-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 20) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 21) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               XY-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 22) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 23) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               XZ-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 24) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 25) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               YZ-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 26) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
            write(6,*) ' '
            write(6,99)
            write(6,99)
            write(6,*) ' '
         endif

      endif

 97   format(5x,"Oscillator Strength:",3x,F12.8)
 98   format(23x,A27)
 99   format(2x,66("-"))
100   format(5x,"Excitation Energy:",5x,F12.8," a.u.",5x,F12.8," eV")
101   format(5x, A5," Transition Moment:",23x,F12.8)
102   format(5x,"Dipole Strength:",31x,F12.8)
103   format(5x,"Oscillator Strengh:",28x,F12.8)
104   format(5x,A35,12x,F12.8)

      return
      end

