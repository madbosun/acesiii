
c These routines localize the recasting of integer types for the system's iargc
c and getarg routines, which handle the command line arguments for Fortran
c programs.

c#define NO_GETARG







      integer function f_iargc()
      integer*4 iargc
      f_iargc = iargc()
      return
      end

      subroutine f_getarg(ndx,sz)
      integer ndx, ierr
      character*(*) sz
      integer*4 i, i2, i3
      i = ndx
      call getarg(i,sz)
      return
      end


