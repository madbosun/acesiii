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
      subroutine get_restart_status(x, nindex, type, bval,
     *                              eval, bdim, edim)
c---------------------------------------------------------------------------
c   Returns 1 or 0 in the scalar argument, depending on whether the current
c   SIAL program has been restarted or not.  If the SIAL program has been
c   restarted, a 1 is returned only for the first call to get_restart_status,
c   all other calls to this instruction will return a 0.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'checkpoint_data.h'

      double precision x
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)

      double precision val

      if (restart_status) then
         val = 1.0
         restart_status = .false.
      else
         val = 0.0
      endif

      x = val
      return
      end
