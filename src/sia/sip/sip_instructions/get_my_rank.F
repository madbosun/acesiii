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
      subroutine get_my_rank(x, nindex, type, bval,
     *                              eval, bdim, edim)
c---------------------------------------------------------------------------
c   Returns the rank of the current processor in the scalar variable.
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

      x = me   ! save the rank in the scalar table.
      return
      end
