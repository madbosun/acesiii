      subroutine gridprmval (
     &     natm, nprm, maxshl, 
     &     nshlatm, nprmshl, angfct,
     &     atmcrd, alpha, pntcrd, prmval,IGRID,NGRID)
c
c The subroutine evaluates the values of the AO primitives at
c a given grid point.
c
c Igor Schweigert, Oct 2002
c $Id: gridprmval.F,v 1.3 2005/04/21 15:42:59 igor Exp $
c 
      implicit none
c
c     Arguments
c
      integer
     &     natm, nprm, maxshl, 
     &     nshlatm (natm), nprmshl (maxshl, natm), angfct (3, maxshl)
c
      double precision
     &     atmcrd (3, natm), alpha (nprm), pntcrd (3), 
     &     prmval (NGRID,nprm)
c
c     Local variables
c
      integer
     &     iatm, ishl, iprm, i,IGRID,NGRID
c
      double precision
     &     dist (3), dst, pre, dexp
c
c     Calculate the values of the primitives at the given point in
c     space: <> Loop over the atoms (the centers of primitives). <>
c     Calculate the distance from the current atom to the point. <> Loop
c     over angular momentum shells on the current atom. <> Calculate the
c     angular momentum prefactor. <> Loop over the primitives of the
c     current shell and atom and calculate their values.
c
      iprm = 1
      do iatm = 1, natm
c     
         dist (1) = pntcrd (1) - atmcrd (1, iatm)
         dist (2) = pntcrd (2) - atmcrd (2, iatm)
         dist (3) = pntcrd (3) - atmcrd (3, iatm)
         dst = dist (1)**2.d0 + dist (2)**2.d0 + dist (3)**2.d0
c
         do ishl = 1, nshlatm (iatm)
c   

            pre =
     &           dist (1)**angfct (1, ishl)*
     &           dist (2)**angfct (2, ishl)*
     &           dist (3)**angfct (3, ishl)
            
c
            do i = 1, nprmshl (ishl, iatm)
          prmval (IGRID,iprm) = pre * dexp (-1.d0 * alpha (iprm) * dst)
               iprm = iprm + 1
            enddo
         enddo
c     
      enddo
c     
      return
      end
