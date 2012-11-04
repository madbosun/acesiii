
c This routine expands a packed super-matrix dA(p<q,*) to dB(pq,*).
c NOTE: dB may alias dA in which case the matrix is expanded "in place".

c CODED JG JUNE/90
c rewritten Q1 2002 ADY

c INPUT
c int irp_pq   : irrep of dA bra (pq pairs)
c int pop(*)   : p/q population vector (orbitals per irrep)
c int max_pq   : distribution size of dB (number of rows)
c int max_pltq : distribution size of dA (number of rows)
c int nCols    : distributions in dA and dB (number of columns)
c double dA(max_pltq,nCols) : source array

c OUTPUT
c double dB(max_pq,nCols) : destination array (may alias dA)

      subroutine symexp2(irp_pq,pop,max_pq,max_pltq,nCols,dB,dA)
      implicit none

c ARGUMENTS
      integer irp_pq, pop(*), max_pq, max_pltq, nCols
      double precision dA(max_pltq,nCols), dB(max_pq,nCols)

c INTERNAL VARIABLES
      integer pq
      integer p, irp_p, max_p
      integer q, irp_q, max_q
      integer off_pltq(8), off_pq(8)
      integer iCol, ndx1, ndx2
      integer iTmp

c COMMON BLOCKS
#include "syminf.com" /* for nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert irp_pq is properly bound
      if ((irp_pq.lt.1).or.(nirrep.lt.irp_pq)) then
         print *, '@SYMEXP2: Assertion failed.'
         print *, '   irp_pq = ',irp_pq
         iTmp = 1
      end if
c   o assert array dimensions are whole
      if ((max_pq.lt.0).or.(max_pltq.lt.0).or.(nCols.lt.0)) then
         print *, '@SYMEXP2: Assertion failed.'
         print *, '   max_pltq = ',max_pltq
         print *, '   max_pq   = ',max_pq
         print *, '   nCols    = ',nCols
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif

      if ((max_pq.lt.0).or.(nCols.lt.0)) return

c ----------------------------------------------------------------------

c   o expand dA into dB (backwards)
      do iCol = nCols, 1, -1
      do pq = max_pltq, 1, -1
         dB(pq,iCol) = dA(pq,iCol)
      end do
      end do

      if (irp_pq.eq.1) then

c      o fill in tables of offsets
         off_pltq(1) = 0
         off_pq(1)   = 0
         do irp_q = 1, nirrep-1
            max_q = pop(irp_q)
            off_pltq(irp_q+1) = off_pltq(irp_q) + max_q*(max_q-1)/2
            off_pq(irp_q+1)   = off_pq(irp_q)   + max_q*max_q
         end do

         do iCol = 1, nCols

            do irp_q = nirrep, 1, -1
               max_q = pop(irp_q)
            if (max_q.gt.1) then

c            o expand the upper-triangle
                  ndx1 = 1 + off_pltq(irp_q) + (max_q-1)*(max_q-2)/2
                  ndx2 = 1 + off_pq(irp_q)   + max_q*(max_q-1)
               do q = max_q-1, 1, -1
               do p = q-1, 0, -1
                  dB(ndx2+p,iCol) = dB(ndx1+p,iCol)
               end do
                  ndx1 = ndx1 - q + 1
                  ndx2 = ndx2 - max_q
               end do

c            o zero the diagonal
                  iTmp = ndx2
               do q = 1, max_q
                  dB(iTmp,iCol) = 0.d0
                  iTmp = iTmp + max_q + 1
               end do

c            o fill in the lower-triangle
               do q = 1, max_q-1
                  ndx1 = ndx2 + max_q + q - 1
               do p = q, max_q-1
                  dB(ndx2+p,iCol) = -dB(ndx1,iCol)
                  ndx1 = ndx1 + max_q
               end do
                  ndx2 = ndx2 + max_q
               end do

c           else if (max_q.lt.2) then
            else

               if (max_q.eq.1) dB(1+off_pq(irp_q),iCol) = 0.d0

c           end if (max_q.gt.1) then
            end if
c           end do irp_q = nirrep, 1, -1
            end do

c        end do iCol = 1, nCols
         end do

c     else if (irp_pq.ne.1) then
      else

c      o fill in tables of offsets
         off_pltq(1) = 0
         off_pq(1)   = 0
         do irp_q = 1, nirrep-1
            irp_p = dirprd(irp_q,irp_pq)
            if (irp_p.lt.irp_q) then
               max_p = pop(irp_p)
               max_q = pop(irp_q)
               iTmp = max_p*max_q
               off_pltq(irp_q+1) = off_pltq(irp_q) + iTmp
               off_pq(irp_q+1)   = off_pq(irp_q)   + iTmp
            else
               max_q = pop(irp_q)
               max_p = pop(irp_p)
               off_pltq(irp_q+1) = off_pltq(irp_q)
               off_pq(irp_q+1)   = off_pq(irp_q) + max_p*max_q
            end if
         end do

         do iCol = 1, nCols

c         o expand the p<q elements
            do irp_q = nirrep, 1, -1
               irp_p = dirprd(irp_q,irp_pq)
               if (irp_p.lt.irp_q) then
                  max_p = pop(irp_p)
                  max_q = pop(irp_q)
                  ndx1 = 1 + off_pltq(irp_q)
                  ndx2 = 1 + off_pq(irp_q)
                  do pq = (max_p*max_q)-1, 0, -1
                     dB(ndx2+pq,iCol) = dB(ndx1+pq,iCol)
                  end do
               end if
            end do

c         o fill in the q<p elements
            do irp_q = 1, nirrep
               irp_p = dirprd(irp_q,irp_pq)
               if (irp_q.lt.irp_p) then
                  max_q = pop(irp_q)
                  max_p = pop(irp_p)
                     ndx2 = 1 + off_pq(irp_q)
                     iTmp = 1 + off_pq(irp_p)
                  do q = 1, max_q
                     ndx1 = iTmp
                  do p = 1, max_p
                     dB(ndx2,iCol) = -dB(ndx1,iCol)
                     ndx1 = ndx1 + max_q
                     ndx2 = ndx2 + 1
                  end do
                     iTmp = iTmp + 1
                  end do
               end if
            end do

c        end do iCol = 1, nCols
         end do

c     end if (irp_pq.eq.1)
      end if

      return
c     end subroutine symexp2
      end

