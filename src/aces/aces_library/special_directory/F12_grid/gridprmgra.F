      subroutine gridprmgra (
     &     natm, nprm, maxshl, 
     &     nshlatm, nprmshl, angfct,
     &     atmcrd, alpha, pntcrd, prmgra,IGRID,NGRID)
c
c The subroutine evaluates the values of primitive gradients at
c a given grid point.
c
c Igor Schweigert, Dec 2002
c $Id: gridprmgra.F,v 1.2 2005/05/03 00:54:03 igor Exp $
c 
      implicit none
c
c     Arguments
c
      integer
     &     natm, nprm, maxshl,IGRID,NGRID, 
     &     nshlatm (natm), nprmshl (maxshl, natm), angfct (3, maxshl)
c
      double precision
     &     atmcrd (3, natm), alpha (nprm), pntcrd (3), 
     &     prmgra (NGRID,nprm,3)
c
c     Local variables
c
      integer
     &     iprm, iatm, ishl, ixyz, i, pow0 (3,3), pow1 (3,3)
c
      double precision
     &     dist (3), dst, e, pre0 (3), pre1 (3)
c
      data 
     &     pow0 /-1,  0,  0,  0, -1,  0,  0,  0, -1/,
     &     pow1 / 1,  0,  0,  0,  1,  0,  0,  0,  1/
c
      iprm = 1
      do iatm = 1, natm
c
c     Calculate the gradients of the primitives: <> Loop over atoms. <>
c     Calculate the distance from the current atom to the grid point. <>
c     Loop over primitive shells on the atom. <> Calculate the
c     exponent-independent prefactors. <> Loop over primitives in the
c     current shell and calculate the gradients.
c
         dist (1) = pntcrd (1) - atmcrd (1, iatm)
         dist (2) = pntcrd (2) - atmcrd (2, iatm)
         dist (3) = pntcrd (3) - atmcrd (3, iatm)
         dst = dist (1)**2 + dist (2)**2 + dist (3)**2
c
         do ishl = 1, nshlatm (iatm)
            do ixyz = 1, 3
               pre0 (ixyz) = 0d0
               if (angfct (ixyz,ishl) .ne. 0 ) pre0 (ixyz) = 
     &              dble(angfct (ixyz,ishl)) *
     &              dist (1)**(angfct(1,ishl)+pow0(1,ixyz))*
     &              dist (2)**(angfct(2,ishl)+pow0(2,ixyz))*
     &              dist (3)**(angfct(3,ishl)+pow0(3,ixyz))
               pre1 (ixyz) = -2d0 *       
     &              dist (1)**(angfct(1,ishl)+pow1(1,ixyz))*
     &              dist (2)**(angfct(2,ishl)+pow1(2,ixyz))*
     &              dist (3)**(angfct(3,ishl)+pow1(3,ixyz))
            enddo
c
            do i = 1, nprmshl (ishl, iatm)
               e = dexp ( - alpha (iprm) * dst)
               do ixyz = 1, 3
                  prmgra (IGRID,iprm,ixyz) = e * (
     &                 pre0 (ixyz) + pre1 (ixyz) * alpha (iprm))
               enddo
               iprm = iprm + 1
            enddo
c     
         enddo
      enddo
c     
      return
      end
