
c This routine initializes the syminf common block.
c The code for the direct product table was taken straight from VMOL.

      subroutine aces_com_syminf
      implicit none

c INTERNAL VARIABLES
      integer i, j, k
      integer zMin, zMax
      logical bTmp

c COMMON BLOCKS
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

c INLINE FUNCTIONS
      integer ka
      ka(i,k) = (i-1)/(2**(k-1)) - 2*((i-1)/(2**k))

c ----------------------------------------------------------------------

c   o nirrep
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)

c   o irrep a/b
c???      norbs=nocco(1)+nvrto(1)
c???      call getrec(-1,'JOBARC','IRREPALP',norbs,irrepa)
c???      call getrec(-1,'JOBARC','IRREPBET',norbs,irrepb)

c   o dirprd
      do i = 1, 8
      do j = 1, 8
         zMin = 1
         zMax = 1
         do k = 1, 3
            zMin = zMin + ( min(ka(i,k),ka(j,k)) * 2**(k-1) )
            zMax = zMax + ( max(ka(i,k),ka(j,k)) * 2**(k-1) )
         end do
         dirprd(i,j) = 1 + zMax - zMin
      end do
      end do


      return
c     end subroutine aces_com_syminf
      end

