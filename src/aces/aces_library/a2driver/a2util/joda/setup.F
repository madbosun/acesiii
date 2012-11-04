      SUBROUTINE SETUP(HES, GRD, HESMOD, GRDMOD, IQFIX, STPMAX, LUOUT)
C
C Build the Hessian and Gradient vectors symmetry coordintes (
C only those internal that are optimized is taking in to account)
C Also do all other initializations.
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "coord.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT, INR, IVEC, IDIE, ICURVY, IMXSTP, ISTCRT, 
     &                IVIB, ICONTL, IRECAL, INTTYP, IDISFD, IGRDFD, 
     &                ICNTYP, ISYM, IBASIS, XYZTol

C
      DIMENSION HES(NXM6, NXM6), GRD(NXM6), HESMOD(NOPT, NOPT),
     &          GRDMOD(NOPT), IQFIX(3*MXATMS, 3*MXATMS)
C
C Let's update the number of cycles right away!
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The number of opt. cycles:", NCYCLE 
#endif
C
      NCYCLE = NCYCLE + 1 
C
      DO 311 I=1, NUNIQUE
         DO 312 J=1, NEQ(IUNIQUE(I))
            IQFIX(IUNIQUE(I), J)=IEQUIV(I,J)
 312     CONTINUE
 311  CONTINUE
C  
      DO 19 I=1, NOPT
         DO 20 J=1, NOPT
C
            HESMOD(I,J)=HES(NOPTI(I), NOPTI(J))
C
            DO 29 K=1, NEQ(NOPTI(I))
               HESMOD(I,J)=HESMOD(I,J)+HES(IQFIX(NOPTI(I),K),NOPTI(J))
 29         CONTINUE
C
            DO 39 L=1, NEQ(NOPTI(J))
               HESMOD(I,J)=HESMOD(I,J)+HES(NOPTI(I),IQFIX(NOPTI(J),L))
 39         CONTINUE
C
            DO 49 K=1, NEQ(NOPTI(I))
               DO 50 L=1, NEQ(NOPTI(J))
                  HESMOD(I, J)=HESMOD(I, J) + HES(IQFIX(NOPTI(I), K),
     &                        IQFIX(NOPTI(J), L))
 50            CONTINUE
 49         CONTINUE
C
            ZD=DSQRT(DFLOAT((NEQ(NOPTI(I))+1)*(NEQ(NOPTI(J))+1)))
            HESMOD(I,J)=HESMOD(I,J)/ZD
C
 20      CONTINUE
 19   CONTINUE
C
#ifdef _DEBUG_LVLM1
          Write(6,*)
          WRITE(LUOUT,115)
 115      FORMAT(T3,' Hessian matrix in totally symmetric symmetry ',
     &      'coordinates: ')
          IF(ICURVY.EQ.1) WRITE(LuOut,117)
 117      FORMAT(T3,'            (Curvilinear representation) ')
          CALL OUTPUT(HESMOD, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
#endif
C     
C Form modified Gradient vector and symmetrize it. 
C
      DO 1177 I=1, NOPT
         GRDMOD(I)=DSQRT(DFLOAT(NEQ(NOPTI(I))+1))*GRD(NOPTI(I))
 1177 CONTINUE
C
#ifdef _DEBUG_LVLM1
      Write(6,*) 
      Write(6,*)  "The modified gradients"
      Write(6, "(3F10.5)") (GRDMOD(I), I=1, NOPT)
      Write(6,*)
#endif
C
C The maximum step size can be controlled by the user by setting
C the flag MAX_STEP. The default value for maximum is set 300.
C IMXSTP is first set in mkvmol.F (strange!). We need to move
C these things into a one place that is visible to the developers.
      IF (IMXSTP .EQ. 0) IMXSTP = 300
      STPMAX = DBLE(IMXSTP)/1000.D0
C
      RETURN
      END




