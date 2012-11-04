
c This routine antisymmetrizes an unpacked super-matrix dA by
c
c    dA(pq,*) = dA(pq,*) - dA(qp,*)
c
c provided p and q are from the same population vector. The output
c array is not packed p<q.

c CODED SEPTEMBER/90 JG
c rewritten Q1 2002 ADY

c INPUT
c int irp_pq : irrep of dA bra (pq pairs)
c int pop(*) : p/q population vector (orbitals per irrep)
c int max_pq : distribution size of dA (number of rows)
c int nCols  : distributions in dA (number of columns)
c double dScr(*) : [OBSOLETE]
c int iScr(*)    : [OBSOLETE]

c INPUT/OUTPUT
c double dA(max_pq,nCols) : array to be antisymmetrized

      subroutine assym2a(irp_pq,pop,max_pq,nCols,dA,dScr,iScr)
      implicit none

c ARGUMENTS
      integer irp_pq, pop(*), max_pq, nCols, iScr(*)
      double precision dA(max_pq,nCols), dScr(*)

c INTERNAL VARIABLES
      integer off_pq(8)
      integer p, irp_p, max_p
      integer q, irp_q, max_q
      logical do_irp_q
      integer iCol, ndx, ndx_up, ndx_lo
      integer iTmp
      double precision dTmp

c COMMON BLOCKS
#include "syminf.com" /* for nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert irp_pq is properly bound
      if ((irp_pq.lt.1).or.(nirrep.lt.irp_pq)) then
         print *, '@ASSYM2A: Assertion failed.'
         print *, '   irp_pq = ',irp_pq
         iTmp = 1
      end if
c   o assert rows and columns are whole
      if ((max_pq.lt.0).or.(nCols.lt.0)) then
         print *, '@ASSYM2A: Assertion failed.'
         print *, '   max_pq = ',max_pq
         print *, '   nCols  = ',nCols
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif

      if ((max_pq.lt.0).or.(nCols.lt.0)) return

c ----------------------------------------------------------------------

      if (irp_pq.eq.1) then

         do iCol = 1, nCols

            ndx = 1

            do irp_q = 1, nirrep
               max_q = pop(irp_q)
            if (max_q.gt.1) then

c            o antisymmetrize the off-diagonal elements
                  ndx_up = ndx + max_q
               do q = 1, max_q-1
                  ndx_lo = ndx + q
               do p = 0, q-1
                  dTmp = dA(ndx_up+p,iCol) - dA(ndx_lo,iCol)
                  dA(ndx_up+p,iCol) =  dTmp
                  dA(ndx_lo,iCol)   = -dTmp
                  ndx_lo = ndx_lo + max_q
               end do
                  ndx_up = ndx_up + max_q
               end do

c            o zero the diagonal elements (while incrementing ndx)
               do q = 1, max_q
                  dA(ndx,iCol) = 0.d0
                  ndx = ndx + max_q + 1
               end do
               ndx = ndx - max_q

c           else if (max_q.lt.2) then
            else

               if (max_q.eq.1) then
                  dA(ndx,iCol) = 0.d0
                  ndx = ndx + 1
               end if

c           end if (max_q.gt.1)
            end if
c           end do irp_q = 1, nirrep
            end do

c        end do iCol = 1, nCols
         end do

c     else if (irp_pq.ne.1) then
      else

         off_pq(1) = 0
         do irp_q = 1, nirrep-1
            irp_p = dirprd(irp_q,irp_pq)
            off_pq(irp_q+1) = off_pq(irp_q) + pop(irp_p)*pop(irp_q)
         end do

         do iCol = 1, nCols

            do irp_q = 1, nirrep
               irp_p = dirprd(irp_q,irp_pq)
               if (irp_p.lt.irp_q) then
                  max_p = pop(irp_p)
                  max_q = pop(irp_q)
                  do_irp_q = ((max_p.ne.0).and.(max_q.ne.0))
               else
                  do_irp_q = .false.
               end if
            if (do_irp_q) then

                  ndx_lo = 1 + off_pq(irp_p)
                  ndx    = 1 + off_pq(irp_q)
               do q = 1, max_q
                  ndx_up = ndx_lo
               do p = 1, max_p
                  dTmp = dA(ndx,iCol) - dA(ndx_up,iCol)
                  dA(ndx,iCol)    =  dTmp
                  dA(ndx_up,iCol) = -dTmp
                  ndx_up = ndx_up + max_q
                  ndx    = ndx    + 1
               end do
                  ndx_lo = ndx_lo + 1
               end do

c           end if (do_irp_q)
            end if
c           end do irp_q = 1, nirrep
            end do

c        end do iCol = 1, nCols
         end do

c     end if (irp_pq.eq.1)
      end if

      return
c     end subroutine assym2a
      end

