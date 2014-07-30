      SUBROUTINE TRANS_SPH2CART(DMAT_IN_OUT, TRANSFRM, TMP, NBAS, NBASP)
C
C This little subroutine makes a in-place transforamtion of a 
C matrix in spherical harmonic basis to a Cartesian basis. 
C
C   DMAT_IN_OUT : Matrix that need to be transformed. Note that the 
C                transformed matrix is also returned in MAT_IN_OUT
C   
C   TRANSFRM    : The transformation matrix is passed in as a input
C                argument
C
C   TMP         : A scratch array to keep the intermediates
C
C  NBASP, NBAS  : Spherical and Cartesian dimensions NBAS > NBASP
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      
      DIMENSION DMAT_IN_OUT(NBASP*NBASP), TRANSFRM(NBAS*NBASP),
     &          TMP(NBAS*NBAS)
      
      DATA ONE, ZLICH /1.0D00, 0.0D00/
C 
      CALL XGEMM('N', 'N', NBAS, NBASP, NBASP, ONE, TRANSFRM, NBAS, 
     &            DMAT_IN_OUT, NBASP, ZLICH, TMP, NBAS)
      
      CALL XGEMM('N', 'T', NBAS, NBAS, NBASP, ONE, TMP, NBAS,
     &           TRANSFRM, NBAS, ZLICH, DMAT_IN_OUT, NBAS)

      RETURN
      END
