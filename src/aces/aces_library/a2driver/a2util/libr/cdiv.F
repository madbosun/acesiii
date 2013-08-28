
c complex division, (cr,ci) = (ar,ai)/(br,bi)

      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
      double precision s,ars,ais,brs,bis,s_inv
c old
c      s = dabs(br) + dabs(bi)
c      ars = ar/s
c      ais = ai/s
c      brs = br/s
c      bis = bi/s
c      s = brs**2 + bis**2
c      cr = (ars*brs + ais*bis)/s
c      ci = (ais*brs - ars*bis)/s
c new
      s_inv = 1.0d0 / ( dabs(br) + dabs(bi) )
      ars = ar * s_inv
      ais = ai * s_inv
      brs = br * s_inv
      bis = bi * s_inv
      s_inv = 1.0d0 / ( brs*brs + bis*bis )
      cr = (ars*brs + ais*bis) * s_inv
      ci = (ais*brs - ars*bis) * s_inv
c end
      return
      end
