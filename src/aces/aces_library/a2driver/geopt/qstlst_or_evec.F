
      SUBROUTINE QSTLST_OR_EVEC(QSTLST_TANGENT, GRDMOD, HESMOD, DIAGHES,
     &                          HES, SCRATCH, IMODE, QSTLST_CLIMB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL QSTLST_CLIMB
      PARAMETER (THRESHOLD = 5.0D-2, WGHT= 1.0D-3)

#include "mxatms.par"
#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DATA HALF, ONE /0.50D0, 1.0D0 /

      DIMENSION HESMOD(NOPT, NOPT), GRDMOD(NOPT), QSTLST_TANGENT(NOPT),
     &          SCRATCH(NX*NX), DIAGHES(NOPT, NOPT),
     &          HES(NOPT, NOPT)

      CALL ZERO(SCRATCH, NOPT)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") "The QST/LST tangent in QSTLST_OR_EVEC"
      Write(6,"(3(1x,F12.6))") (QSTLST_TANGENT(I), I= 1, NOPT)
#endif
      CALL XGEMM('N','N', NOPT, 1, NOPT, 1.0D0, HES, NOPT,
     &            QSTLST_TANGENT, NOPT, 0.0D0, SCRATCH, NOPT)

      EIGVALUE = XDDOT(NOPT, QSTLST_TANGENT, 1, SCRATCH, 1)
      GRDVALUE = XDDOT(NOPT, QSTLST_TANGENT, 1, GRDMOD,  1)

      EST_STEP = (EIGVALUE + DSQRT(EIGVALUE**2 + 4.0D0*(GRDVALUE)
     &            **2))*HALF
      EST_STEP =  -GRDVALUE/(EST_STEP - EIGVALUE)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)")"Relevent parameter that is used to decide a climb"
      Write(6,"(a,2(1x,F12.6))") "The EIGVALUE,GRDVALUE : ", EIGVALUE, 
     &                            GRDVALUE
      Write(6,"(a,F12.6)") "The estimated step size     :",EST_STEP
      Write(6,*)
      Write(6,"(a)") "The HT vector"
      Write(6,"(3(1x,F12.6))") (Scratch(i), i=1, NOPT)
#endif 
C
C Let's check that the estimated step size is greater than 0.05 au.
C If so, then follow the QST or LST tangent (I am interpreting Schlegel's
C "estimated displacement along the tangent vector is greater than the 0.05".
C At this point I am not sure whether he meant the largest absolute
C displacement. Also cap the QST/LST to first four cycle (also rec. by Schlegel)
C
      QSTLST_CLIMB = ((DABS(EST_STEP).GE.THRESHOLD.AND.NCYCLE.LE.4).OR.
     &                NCYCLE.LE.2)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,l)") "The climbing Phase if true", QSTLST_CLIMB
#endif
C
C Compute the gradient along the tangent vector, copy the tangent
C vector to the eigenvector matrix (to the first eigenvector) and the
C eigenvalue (EIGVALUE) to the first eigenvalue and set the IMODE
C to 1.
C
      IF (QSTLST_CLIMB) THEN
         WEIGHT = ONE - WGHT
         CALL MODFY_HESSIAN(DIAGHES, HESMOD, HES, QSTLST_TANGENT,
     &                      SCRATCH, EIGVALUE, WEIGHT, NOPT)
      END IF
C
C This is just to be safe in subsequent steps.
C
      CALL ZERO(SCRATCH,NOPT*NOPT)

      RETURN
      END

