
C VERY IMPORTANT ROUTINE. INITIALIZES ENVIRONMENT FOR A CRAPS MODULE.
C THIS MUST BE CALLED IN THE FIRST EXECUTABLE STATEMENT IN ALL PROGRAMS.
C
C TO USE THE AUTOMATIC MEMORY ALLOCATION, IENTRY MUST BE SET TO A
C VALUE OTHER THAN -1.

c#define _DEBUG_CRAPSI

      SUBROUTINE CRAPSI(ICORE,IUHF,IENTRY)
      IMPLICIT INTEGER (A-Z)

      INTEGER ICORE(1), tmp
      DOUBLE PRECISION X

c PARAMETERS
      INTEGER IONE
      PARAMETER (IONE=1)

c COMMON BLOCKS
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "crapsx_time.com"
      integer iflag
      common /iflag/ iflag
      data    iflag /0/

c ----------------------------------------------------------------------

#ifdef _TRAPFPE
cYAU Enable floating-point exception handling if the compiler cannot.
      call trapfpe
#endif /* _TRAPFPE */

cjp setup stdout buffering in order to see output immediatelly
      call bufferlines

      NFLAGS  = 100
      NFLAGS2 = 500
      CALL MCHPRM

      call dumpja('I')
c      CALL GETREC(20,'JOBARC','IENDSTAT',IONE,ISTAT)
      CALL CHECKIN
      call wallclock(year_i,mon_i,mday_i,hour_i,min_i,sec_i)

      IF (IENTRY.NE.-1) THEN
         CALL GETREC(20,'JOBARC','IFLAGS  ',NFLAGS,IFLAGS)
         CALL GETREC(20,'JOBARC','IFLAGS2 ',NFLAGS2,IFLAGS2)

C ABACUS CAUSES SOME MEMORY PROBLEMS WHEN REQUIRING THE FULL CORE
C WE THEREFORE REDUCE THE REQUESTED CORE SIZE BY A CERTAIN AMOUNT
C BUT WE MAKE SURE THAT A MINIMUM IS ALWAYS ALLOCATED
C CALL CRAPSI THERFORE WITH IENTRY = 100
         IF (IENTRY.EQ.100) THEN
            IENTRY=0
            MEM=IFLAGS(36)
C         o set to run efficiently in the 8MW queue of a YMP
            MEMMIN=MIN(3500000,IFLAGS(36))
            MEMMAX=MAX(MEMMIN,IFLAGS(36)-3000000)
         ELSE
            MEMMAX=IFLAGS(36)
         END IF

c      o allocate memory once (iflag is only used by BWCC currently)
         if (iflag.eq.0) then
            call aces_malloc(memmax,i0)
            if (i0.eq.-1) then
               write(*,*) '@CRAPSI: Request for ',MEMMAX,
     &                    ' words of memory failed.'
               call errex
            end if
#ifdef _DEBUG_CRAPSI
            write(*,*) '@CRAPSI: memmax = ',memmax
            write(*,*) '         i0     = ',i0
#endif /* _DEBUG_CRAPSI */
            iflag=1
         end if

C Make sure that I0 starts on a double precision word boundary

#ifdef _UNICOS
         I0 = I0 + IINTFP - MOD(I0,IINTFP) - 1
#else
C Coded by Marcel Nooijen and incorporated by Ajith
         IPASS = 0
 100     tmp=(LOC(ICORE(I0)) - LOC(X))
         IF (MOD(tmp, IFLTLN) .NE. 0) THEN
            IF (IPASS .GE. IINTFP) THEN
               WRITE(6, *) '@CRAPSI: Memory alignment failed.'
               CALL ERREX
            END IF
            I0     = I0 + 1
            MEMMAX = MEMMAX - 1
            IPASS  = IPASS + 1
            GO TO 100
         END IF
         tmp = (LOC(ICORE(I0)) - LOC(X))
         IF (MOD(tmp,IFLTLN) .NE. 0) THEN
            WRITE(6,*) '@CRAPSI: Memory alignment problem persists.'
            CALL ERREX
         END IF
#endif /* _UNICOS */

         CALL POPCORN(MEMMAX)
      ELSE
         I0=1
      END IF

      CALL GETSTF(ICORE(I0),IUHF,IENTRY)

      RETURN
      END

