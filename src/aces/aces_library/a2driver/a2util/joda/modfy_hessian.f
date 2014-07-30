
      SUBROUTINE MODFY_HESSIAN(DIAGHES, HESMOD, HES, QSTLST_TANGENT,
     &                         SCRATCH, EIGVALUE, WEIGHT, NOPT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION HESMOD(NOPT, NOPT), QSTLST_TANGENT(NOPT),
     &          SCRATCH(NOPT*NOPT), DIAGHES(NOPT, NOPT),
     &          HES(NOPT, NOPT)
      DATA ONE /1.0D0/, TWO /2.0D0/
C
C Note that at the moment SCRATCH array contains the HxT. Let's label
C that product as A vector. As we can see the A vector is an arbitrary
C direction in the space of the eigenvectors of H.
C
      SCALED_WEIGHT = TWO*EIGVALUE*WEIGHT
      CALL XGEMM('T','N', NOPT, NOPT, 1, -WEIGHT, QSTLST_TANGENT,
     &            NOPT, SCRATCH, NOPT, ONE, HES, NOPT)
      CALL XGEMM('N','T', NOPT, NOPT, 1, -WEIGHT, QSTLST_TANGENT,
     &            NOPT, SCRATCH, NOPT, ONE, HES, NOPT)
      CALL XGEMM('T','N', NOPT, NOPT, 1, SCALED_WEIGHT,
     &            QSTLST_TANGENT, NOPT, QSTLST_TANGENT,
     &            NOPT, ONE, HES, NOPT)
C
C Diagonalize the modified Hessian and print the eigenvalues and
C vectors. These new vectors and the values dictate the climbing
C phase of the search.
C
      CALL DCOPY(NOPT*NOPT, HES, 1, HESMOD, 1)
      CALL EIG(HESMOD, DIAGHES, NOPT, NOPT, 1)

c      Write(6,*) "The new Eigen vectors"
c      Call output(DIAGHES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
c      Write(6,*) "The new Eigen values"
c      Call output(HESMOD, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)

      RETURN
      END

