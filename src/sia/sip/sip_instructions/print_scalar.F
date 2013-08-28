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
      subroutine print_scalar(x, nindex, type, bval,
     *                              eval, bdim, edim)
c---------------------------------------------------------------------------
c   Prints the value of a scalar variable.  The scalar to be printed
c   is defined in the c_result_array field of the op argument.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      double precision x
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      if (me .eq. 0) 
     *      print *,'Task ',me,' Scalar value = ',
     *      x,' at line number ',current_line
      call c_flush_stdout()
      return
      end
