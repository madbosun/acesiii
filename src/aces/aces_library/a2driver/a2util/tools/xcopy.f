















c#define _ABORT_ON_OVERLAP

      subroutine xcopy(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,n
      INTEGER c_adr, z
      external c_adr
      if(n.le.0)return
      if ((incx.eq.1).and.(incy.eq.1)) then
         if (c_adr(dy).eq.c_adr(dx)) return
         call c_memmove(dy,dx,n*8)
      else
         ix = 1
         iy = 1
         if(incx.lt.0)ix = (-n+1)*incx + 1
         if(incy.lt.0)iy = (-n+1)*incy + 1
         do i = 1,n
           dy(iy) = dx(ix)
           ix = ix + incx
           iy = iy + incy
         end do
      end if
      return
      end

