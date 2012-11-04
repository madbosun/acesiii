
#ifdef __fortran9x

! This routine is used to allocate memory.

! iwords : the number of *integer* words requested
! ioff   : the offset (in integer words) relative to icore(1)
!          OR -1 if the allocation failed

SUBROUTINE aces_malloc(iwords,ioff)
   IMPLICIT NONE

   INTEGER, ALLOCATABLE, SAVE :: heap(:)
   INTEGER, INTENT(IN)    :: iwords
   INTEGER, INTENT(INOUT) :: ioff

#include "icore.com"
#include "machsp.com"

   ALLOCATE(heap(iwords),STAT=ioff)
   IF (ioff /= 0) THEN
      ioff = -1 ! leave error handling to the caller
      RETURN
   END IF
   ioff = LOC(heap) - LOC(icore) ! the byte offset
   ioff = 1 + (ioff/iintln)      ! scale the offset to integers
#ifdef _DEBUG
   write(*,*) '@ACES_MALLOC: Allocated ~',iwords*BIT_SIZE(ioff)/8388608,' MB'
#endif

 ! ADY 20000203: The rest is supposedly for Crays, but I can't confirm
 ! this since switching to FORTRAN 95.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! CALL timef(telap) !declare telap as double precision
 ! CALL hpalloc(index,iwords,ierr1,0)
 ! CALL timef(telap)
 ! WRITE(*,*) '@ACES_MALLOC: Allocation of ',iwords,' words of core required ',&
 !  & TELAP/1000.0,' sec.'
 ! IF (IERR1 < 0) THEN
 !    ioff=-1
 !    RETURN
 ! END IF
 ! ioff=1+(index-LOC(icore))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   RETURN
END

#else /* defined(__fortran77) */

c This routine is used to allocate memory.

c iwords : the number of *integer* words requested
c ioff   : the offset (in integer words) relative to icore(1)
c          OR -1 if the allocation failed

      SUBROUTINE aces_malloc(iwords,ioff)
      IMPLICIT NONE

      INTEGER iwords, ioff, ibytes

#include "icore.com"
#include "machsp.com"

#ifdef _UNICOS
      INTEGER index, ierr1
      CALL hpalloc(index,iwords,ierr1,0)
      IF (ierr1.LT.0) THEN
         ioff = -1
         RETURN
      END IF
      ioff = 1 + (index-LOC(icore))
#else
#   ifdef _DEC
C Ajith 05/20/96 for DecAlpha machines. Having to do with the
C different integer representations C and Fortan on DecAlpha.
      INTEGER*8 index, ioff, MALLOC, loc
      ibytes = iwords * iintln
      index = MALLOC(ibytes)
#   else
      INTEGER index, c_malloc
      ibytes = iwords * iintln
      index = c_malloc(ibytes)
#   endif
      IF (index.EQ.0) THEN
         ioff = -1
         RETURN
      END IF
      ioff = index - loc(icore)
      ioff = 1 + (ioff/iintln)
#ifdef _DEBUG
      write(*,*) '@ACES_MALLOC: Allocated ~',rshift(ibytes,20),' MB'
#endif
#endif

      RETURN
      END

#endif /* __fortran9x */

