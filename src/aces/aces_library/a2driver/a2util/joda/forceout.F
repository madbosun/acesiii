      SUBROUTINE FORCEOUT(V)
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
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
 
#include "io_units.par"
 
      DIMENSION V(NXM6)
C   
      WRITE(LuOut,80)
 80   FORMAT(T3,' Internal coordinate forces and energy gradients ',
     &'(atomic units): '
     &,/,T12,'R',T20,'dV/dR',T38,'R',T46,'dV/dR',T64,'R',T72,'dV/dR')

      WRITE(LUOUT, *)
      WRITE(LUOUT,81)(VARNAM(ISQUASH(I)),R(I),V(I),I=1,NXM6)
cYAU 81   FORMAT(3(1X,:'[',A,']',1X,F8.5,1X,F10.7))
 81   FORMAT(3(1X,:'[',A,']',1X,F8.5,1X,F8.5))
      WRITE(LUOUT, *)

      RETURN
      END

