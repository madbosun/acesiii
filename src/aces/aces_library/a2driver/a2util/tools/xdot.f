
      double precision function xdot(n,dx,incx,dy,incy)

      double precision dx(*),dy(*),ddot
      integer incx, incy, n




      integer   int_n, int_incx, int_incy


      if (n.le.0) return

      int_incx = incx
      int_incy = incy
      int_n    = n




      xdot = ddot(int_n,dx,int_incx,dy,int_incy)


      return
      end

