
c This transforms a triangular block of a matrix to a square one
      subroutine mat_trans_tri_sqr(tri,sqr,trilen,sqrlen,sqroff)


















c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler


























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      implicit none


c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/


      integer trilen,sqrlen,sqroff
      double precision tri(trilen),sqr(sqrlen,sqrlen)

      integer i,j,itri,len

      callstack_curr='MAT_TRANS_TRI_SQR'
      len=int(sqrt(real(1+8*trilen))-1)/2
      itri=0
      do i=sqroff,sqroff+len-1
        do j=sqroff,i
          itri=itri+1
          sqr(i,j)=tri(itri)
          sqr(j,i)=tri(itri)
        end do
      end do
      return
      end
