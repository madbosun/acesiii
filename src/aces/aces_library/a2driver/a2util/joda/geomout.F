      SUBROUTINE GEOMOUT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
#include "cbchar.com"
#include "coord.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
#include "io_units.par"
C
      WRITE(LuOut,80)
 80   FORMAT('         -----------------------------------------------'
     &,/,3x,'      Cartesian coordinates corresponding to internal '
     &,/,3x,      '              coordinate input (Angstroms) ')
      CALL DUMPCORD(NATOMS,Q,IATNUM)
c     &,/,3x,      '              coordinate input (Bohr) ',/
c     &,'      ------------------------------------------------------ ',/
c     &,'       Atomic                               ',/
c     &,'       Symbol            X             Y             Z ' )
c      WRITE(LuOut,*)
c     $ '       ------------------------------------------------------ '
C      DO 20 I=1,NX-2,3
C      WRITE(LuOut,90)ZSYM((I/3)+1),Q(I),Q(I+1),Q(I+2)
C90    FORMAT(9X,A,T18,3(F12.7,2X))
C  20  CONTINUE
C      WRITE(LuOut,*)
c     $ '       ------------------------------------------------------ '
c      WRITE(LuOut,*)' '
      RETURN
      END
