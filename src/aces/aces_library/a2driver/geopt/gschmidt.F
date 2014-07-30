      SUBROUTINE GSCHMIDT(VORTH,NSIZE,NDIM,TMP,RESID)
C
C THIS PROJECTS OUT ALL PARTS OF AN INPUT VECTOR (VEC)
C WHICH LIE IN THE SPACE SPANNED BY THE ORTHOGONAL BASIS
C VORTH.
C
C   |v'> = |v> - SUM <i|v> |i>
C                 i 
C
C WHERE THE |i> ARE NORMALIZED BASIS VECTORS FOR THE SPACE VORTH
C
C INPUT:
C       VEC : THE VECTOR WHICH IS TO BE ORTHOGONALIZED TO
C             THE EXISTING BASIS.  *IT IS ASSUMED THAT VEC
C             IS NORMALIZED ON INPUT)
C     VORTH : THE BASIS VECTORS FOR THE EXISTING ORTHOGONAL
C             BASIS
C     NSIZE : THE LENGTH OF THE BASIS VECTORS
C      NDIM : THE DIMENSION OF THE ORTHOGONAL SPACE
C       TMP : A SCRATCH VECTOR OF LENGTH NSIZE
C     RESID : THE NORM OF VEC, AFTER ORTHONALIZATION AND
C             BEFORE NORMALIZATION
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VORTH(NSIZE,NDIM),TMP(NSIZE)
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      write(6,"(a)") 'Non-orthogonal vectors in gschmidt'
      call output(vorth,1,nsize,1,ndim,nsize,ndim,1)
      do i = 1, ndim
         TMP(i)=XDDOT(NSIZE,VORTH(1,I),1,VORTH(1,I),1)
      enddo  
      Write(6,*)
      write(6,"(a)") 'Old overlap in gschmidt'
      Write(6, "(6(1x,F12.6))") (TMP(i), i=1, NSIZE)
#endif
C
      DO I=2, NDIM 
         DO J=1, I-1
            FACT=XDDOT(NSIZE,VORTH(1,I),1,VORTH(1,J),1)
            CALL XDAXPY(NSIZE,-FACT, VORTH(1,J), 1,VORTH(1,I),1)
            TMP(I)=XDDOT(NSIZE,VORTH(1,I),1,VORTH(1,I),1)
            CALL XDSCAL(NSIZE, TMP(I), VORTH(1,I), 1)
         ENDDO
      ENDDO 

C
      DO I = 2, NDIM
         TMP(I)=XDDOT(NSIZE,VORTH(1,I),1,VORTH(1,I),1)

         CALL XDSCAL(NSIZE, 1.0D0/DSQRT(TMP(I)), VORTH(1,I), 1)
      ENDDO

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") 'Orthogonal vectors in gschmidt'
      call output(vorth,1,nsize,1,ndim,nsize,ndim,1)
      do i = 1, ndim
         TMP(i)=XDDOT(NSIZE,VORTH(1,I),1,VORTH(1,I),1)
      enddo
      Write(6,*)
      write(6,"(a)") 'Current overlap in gschmidt'
      Write(6, "(6(1x,F12.6))") (Tmp(i), i=1, NSIZE)
#endif
C
      RETURN
      END

