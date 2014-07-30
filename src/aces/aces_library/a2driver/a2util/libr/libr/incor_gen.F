
C This routine loads lists into ICORE(I0), increments I0, and decrements
C ICRSIZ. I0 and ICRSIZ do not have to point to the /ISTART/ common
C block, but they must be relative to ICORE(1) in the blank common
C block. The logic is also reentrant, so it will not load lists that are
C already in memory.

#include "disttype.h"

      SUBROUTINE INCOR_GEN(I0,ICRSIZ,IUHF)
      IMPLICIT NONE
      INTEGER I0,ICRSIZ,IUHF
#include "icore.com"
#include "machsp.com"
#include "flags.com"
#include "sympop.com"
      integer           quikget(10,500)
      common /auxcache/ quikget
      integer iLim, iGrp, iFam, nRows, nCols, nDbls, nInts
      integer  aces_list_rows, aces_list_cols
      external aces_list_rows, aces_list_cols

      IF (IFLAGS(35).EQ.0) RETURN

      IF (IFLAGS(35).EQ.6) THEN
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,1,500,1)
         RETURN
      END IF

      iLim = 0
      IF (IFLAGS(35).EQ.1) THEN
c      o NOABCD
         iLim=max(iLim,irpdpd(1,h_DistType_Ab)*irpdpd(1,h_DistType_Ai))
         iLim=max(iLim,irpdpd(1,h_DistType_Ab)*irpdpd(1,h_DistType_Ia))
         iLim=max(iLim,irpdpd(1,h_DistType_AB)*irpdpd(1,h_DistType_AI))
         iLim=max(iLim,irpdpd(1,h_DistType_ab)*irpdpd(1,h_DistType_ai))
      ELSE IF (IFLAGS(35).EQ.4) THEN
c      o NOABCI
         iLim=max(iLim,irpdpd(1,h_DistType_AI)*irpdpd(1,h_DistType_ai))
         iLim=max(iLim,irpdpd(1,h_DistType_AI)*irpdpd(1,h_DistType_AI))
         iLim=max(iLim,irpdpd(1,h_DistType_ai)*irpdpd(1,h_DistType_ai))
      ELSE IF (IFLAGS(35).EQ.5) THEN
c      o NOABIJ
         iLim=max(iLim,irpdpd(1,h_DistType_Ij)*irpdpd(1,h_DistType_Ia))
         iLim=max(iLim,irpdpd(1,h_DistType_Ij)*irpdpd(1,h_DistType_Ai))
         iLim=max(iLim,irpdpd(1,h_DistType_IJ)*irpdpd(1,h_DistType_IA))
         iLim=max(iLim,irpdpd(1,h_DistType_ij)*irpdpd(1,h_DistType_ia))
      END IF
      IF (iLim.eq.0) RETURN
      iLim = iintfp*iLim

      do iFam = 1, 500
      do iGrp = 1, 10
         if (quikget(iGrp,iFam).eq.0) then
            nCols = aces_list_cols(iGrp,iFam)
            nRows = aces_list_rows(iGrp,iFam)
            nInts = iintfp*nRows*nCols
            if (0.lt.nInts.and.nInts.le.iLim.and.nInts.lt.iCrSiz) then
c               print *, '@INCOR_GEN: Loading ',iGrp,iFam
               call getlst(icore(i0),1,nCols,1,iGrp,iFam)
               quikget(iGrp,iFam) = i0
               i0     = i0     + nInts
               iCrSiz = iCrSiz - nInts
            end if
         end if
      end do
      end do

      RETURN
      END

