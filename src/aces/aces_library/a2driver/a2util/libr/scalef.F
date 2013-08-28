
      subroutine scalef(d,ntot,n)
      implicit double precision(a-h,o-z)
      dimension d(ntot),n(8)
      integer dirprd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)

      ioff = 0
      do irrep = 1, nirrep
         max_n = n(irrep)
         itmp = 1
         do i = 1, max_n
            d(ioff+itmp) = d(ioff+itmp) * 2.0d0
            itmp = itmp + i + 1
         end do
         ioff = ioff + ishft(max_n*(max_n+1),-1)
      end do
      return 
      end

