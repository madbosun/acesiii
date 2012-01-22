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
      subroutine remove_atom_rud2(x, nindex, type, bval,
     *                              eval, bdim, edim)
c--------------------------------------------------------------------------
c Takes a four-dimensional array array(a,b,c,d) and removes the elements
c in which atom(d) = atom(a,b,or c).    
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      if (nindex .eq. 4) then

         call doremove_atom_rud2(x, bval(1),eval(1),type(1),
     *                   bval(2),eval(2),type(2),
     *                   bval(3),eval(3),type(3),
     *                   bval(4),eval(4),type(4))
      else
         print *,'Error in remove_atom_rud1: Result array ',
     *      'must have either 4 indices.'
         print *,'array has ',nindex,' indices.'
         call abort_job()
      endif

      return
      end

      subroutine doremove_atom_rud2(x, a1,a2,atype,
     *                      i1,i2,itype,b1,b2,btype,
     *                      j1,j2, jtype)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h' 

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      integer a,b,c,d,i,j, m1, m2, n1, n2, r1, r2, s1, s2 
      integer m, n, r, s, aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2 
      double precision val

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through to determine
c   the atomic labels .
c-----------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, i1, n1)
      call lookup_shell(end_nfps, nshells, i2, n2)
      call lookup_shell(end_nfps, nshells, b1, r1)
      call lookup_shell(end_nfps, nshells, b2, r2)
      call lookup_shell(end_nfps, nshells, j1, s1)
      call lookup_shell(end_nfps, nshells, j2, s2)

      do m = m1, m2 
            if (m .eq. 1) then
               aa1 = 1
            else
               aa1 = end_nfps(m-1) + 1
            endif
            aa2 = end_nfps(m) 
         do n = n1, n2
            if (n .eq. 1) then
               bb1 = 1
            else
               bb1 = end_nfps(n-1) + 1
            endif
            bb2 = end_nfps(n)
         do r = r1, r2
            if (r .eq. 1) then
               cc1 = 1
            else
               cc1 = end_nfps(r-1) + 1
            endif
            cc2 = end_nfps(r)
         do s = s1, s2
               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)  

             if ((atom(s) .eq. atom(m)) .or. (atom(s) .eq. atom(n)) .or. 
     *           (atom(s) .eq. atom(r))) then 

             do a = aa1, aa2  
             do b = bb1, bb2  
             do c = cc1, cc2  
             do d = dd1, dd2  
                x(a,b,c,d) = 0.0d0 
             enddo 
             enddo 
             enddo 
             enddo 

             endif 

         enddo ! m 
         enddo ! n 
         enddo ! r 
      enddo ! s 


      return
      end


