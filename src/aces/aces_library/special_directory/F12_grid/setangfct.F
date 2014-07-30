      subroutine setangfct (
     &     maxshl, angfct)
c
c This routine calculate the angular factors (lx, ly, lz) of Gaussian
c primitives x**lx * y**ly * z**lz * exp (..) for each shell from 1 to
c MAXSHL. The loop is taken from KS.SRC/INTGRT/INITGRID.
c
c Igor Schweigert, Jul 2004
c $Id: setangfct.F,v 1.1 2004/12/21 17:41:16 igor Exp $
c
      implicit none
c
c     Arguments
c
      integer
     &     maxshl, angfct (3, maxshl)
c
c     Local variables
c
      integer
     &     iang, jang, kang, ishl, x, y, z
c
c     Loop over angular momentum and set the angular factors for each
c     shell.
c     
c     angfct (1, *) = 0  100  211000  3221110000  ...
c     angfct (2, *) = 0  010  010210  0102103210  ...
c     angfct (3, *) = 0  001  001012  0010120123  ...
c
      ishl = 1
      do iang = 1, 100
         if (ishl.gt.maxshl) goto 100
         x = iang - 1
         do jang = 1, iang
            y = jang - 1
            z = 0
            do kang = 1, jang
               angfct (1, ishl) = x
               angfct (2, ishl) = y
               angfct (3, ishl) = z

c               write(6,*) 'X=', angfct (1, ishl), ishl
c               write(6,*) 'Y=', angfct (2, ishl), ishl
c               write(6,*) 'Y=', angfct (3, ishl), ishl

               y = y - 1
               z = z + 1
               ishl = ishl + 1
            enddo
            x = x - 1
         enddo
      enddo
c
 100  continue
      return
      end
