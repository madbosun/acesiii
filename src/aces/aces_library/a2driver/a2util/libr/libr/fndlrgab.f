
C EMULATOR OF CRAY SCILIB ROUTINE. RETURNS THE LARGEST ABSOLUTE VALUE IN A
C VECTOR (X).

c There is no UNICOS man page available for this routine.

      double precision function fndlrgab(x,n)
      integer n,ndx,idamax
      double precision x(n)
      external idamax
      if (n.gt.0) then
         ndx = idamax(n,x,1)
         fndlrgab = abs(x(ndx))
      else
         fndlrgab = 0.d0
      end if
      return
      end

