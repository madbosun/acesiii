
      SUBROUTINE MODFY_HESSIAN(DIAGHES, HESMOD, HES, QSTLST_TANGENT,
     &                         SCRATCH, EIGVALUE, WEIGHT, NOPT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "jodaflags.com"

      DIMENSION HESMOD(NOPT, NOPT), QSTLST_TANGENT(NOPT),
     &          SCRATCH(NOPT), DIAGHES(NOPT, NOPT),
     &          HES(NOPT, NOPT)
      DATA ONE /1.0D0/, TWO /2.0D0/
C
C Note that at the moment SCRATCH array contains the HxT. Let's label
C that product as A vector. As we can see the A vector is an arbitrary
C direction in the space of the eigenvectors of H.
C Impetus for this formula is from Ayala & Schlegel, JCP, 375, 1997.
C Eqn. 7a. 
 
         SCALED_WEIGHT = TWO*EIGVALUE*WEIGHT
        
         CALL XGEMM('T','N', NOPT, NOPT, 1, 1.0D0, QSTLST_TANGENT,
     &               1, SCRATCH, 1, 0.0D0, DIAGHES, NOPT)
         CALL TRANSP(DIAGHES, HESMOD, NOPT, NOPT)
         CALL XDAXPY(NOPT*NOPT, 1.0D0, HESMOD, 1, DIAGHES, 1) 
         CALL XDSCAL(NOPT*NOPT, -WEIGHT, DIAGHES, 1)
         CALL XDAXPY(NOPT*NOPT, 1.0D0, DIAGHES, 1, HES, 1)
         CALL XGEMM('T','N', NOPT, NOPT, 1, SCALED_WEIGHT, 
     &               QSTLST_TANGENT, 1,
     &               QSTLST_TANGENT, 1, 1.0D0, HES, NOPT)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The modified. Hessian"
      Call output(HES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
#endif
C
C Diagonalize the modified Hessian and print the eigenvalues and
C vectors. These new vectors and the values dictate the climbing
C phase of the search.
C
         CALL XDCOPY(NOPT*NOPT, HES, 1, HESMOD, 1)
         CALL EIG(HESMOD, DIAGHES, NOPT, NOPT, 1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The Eigen vectors of modified. Hessian"
      Call output(DIAGHES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
      Write(6,"(a)") "The Eigen values of modified Hessian"
      Write(6, "(3(1x,F12.6)))") (hesmod(i,i), i=1,nopt)
#endif
C
C
C This block of code is obsoleted. 
#ifdef _NO_OBSOLETE
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The Eigen vectors of the orig. Hessian"
      Call output(DIAGHES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
      Write(6,"(a)") "The Eigen values of the orig.  Hessian"
      Write(6, "(6(1x,F12.6))") (hesmod(i,i), i=1,nopt)
#endif
         SCALED_WEIGHT = TWO*EIGVALUE*WEIGHT
         CALL XDCOPY(NOPT, QSTLST_TANGENT, 1, DIAGHES, 1)
         CALL GSCHMIDT(DIAGHES, NOPT, NOPT, SCRATCH, RESID)
C
         CALL XGEMM('T','N', NOPT, NOPT, 1, SCALED_WEIGHT, 
     &               QSTLST_TANGENT, 1,
     &               QSTLST_TANGENT, 1, 1.0D0, HESMOD, NOPT)
      
CSSS      CALL XGEMM('N','N', NOPT, NOPT, NOPT, 1.0D0, DIAGHES, NOPT, 
CSSS     &            HES, NOPT, 0.0D0, SCRATCH, NOPT)
CSSS      CALL XGEMM('N','T', NOPT, NOPT, NOPT, 1.0D0, SCRATCH, NOPT, 
CSSS     &            DIAGHES, NOPT, 0.0D0, HESMOD, NOPT)
      CALL EIG(HESMOD, DIAGHES, NOPT, NOPT, 1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The Eigen vectors of modf. Hessian"
      Call output(DIAGHES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
      Write(6,"(a)") "The Eigen values of modified Hessian"
      Write(6, "(6(1x,F12.6))") (hesmod(i,i), i=1,nopt)
#endif
      Endif
#endif

      RETURN
      END

