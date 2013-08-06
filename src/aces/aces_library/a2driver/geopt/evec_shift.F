      SUBROUTINE EVEC_SHIFT(EIGEN, VEC, SCR, NOPT)
C 
C This routine moves the 3N-6(or 5) translational and rotational
C modes to the end (ie. This eigenvalues and vectors that have
C zero eigenvalues). Ajith Perera, 11/2012 (the date is when
C the comment was added, the routine was written 3 or 4 years 
C prior to that).  

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION EIGEN(NOPT, NOPT), VEC(NOPT, NOPT), SCR(NOPT, NOPT)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Unsorted Eigen values and vectors"
      call output(eigen, 1, nopt, 1, nopt, nopt, nopt, 1)
      call output(vec, 1, nopt, 1, nopt, nopt, nopt, 1)
#endif

      NTR=0
      DO I=NOPT,1,-1
         IF(ABS(EIGEN(I,I)).LT.1.0D-09) NTR=NTR+1
         IF(ABS(EIGEN(I,I)).LT.1.0D-09) ISHIFT=I
      ENDDO
C 
      JSHIFT=ISHIFT-1
      DO J=1,NTR
         JSHIFT=JSHIFT+1
          SCR(8,J)=EIGEN(JSHIFT,JSHIFT)
          DO I=1,NOPT
             SCR(J,I)=VEC(I,JSHIFT)
          ENDDO
      ENDDO
C
      DO J=ISHIFT,NOPT-NTR
        EIGEN(J,J)=EIGEN(J+NTR,J+NTR)
        DO I=1,NOPT
           VEC(I,J)=VEC(I,J+NTR)
        ENDDO
      ENDDO
C
      JTMP=0
      DO J=NOPT-NTR+1,NOPT
         JTMP=JTMP+1
         EIGEN(J,J)=SCR(8,JTMP)
         DO I=1,NOPT
            VEC(I,J)=SCR(JTMP,I)
         ENDDO
      ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Sorted Eigen values and vectors"
      call output(eigen, 1, nopt, 1, nopt, nopt, nopt, 1)
      call output(vec, 1, nopt, 1, nopt, nopt, nopt, 1)
#endif
      
      RETURN
      END


