      SUBROUTINE BUILT_BGMTRX(CARTCOORD, REDUNCO, IREDUNCO,
     &                        TOTREDNCO, TOTNOFBND, TOTNOFANG,
     &                        TOTNOFDIH, NRATMS, 
     &                        BMATRX, GMATRX, EPSILON)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "io_units.par"
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C 
      DIMENSION CARTCOORD(3*NRATMS), IREDUNCO(4, MAXREDUNCO),
     &          BMATRX(TOTREDNCO, 3*NRATMS), 
     &          BTMP(3*MXATMS, 3*MXATMS), 
     &          GMATRX(TOTREDNCO, 3*NRATMS),
     &          EIGVECTORS(TOTREDNCO, TOTREDNCO),
     &          REDUNCO(TOTREDNCO)
         
      DATA MONE /-1.0/
C  
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
      CALL ZERO(BMATRX, TOTREDNCO*(3*NRATMS))
C
      DO 20 IBNDS = 1, TOTNOFBND
C
         IF (IREDUNCO(2, IBNDS) .NE. IZERO) THEN  
C
            CALL  BULT_BNDCRD(CARTCOORD, BMATRX, DISTAB,
     &                        IREDUNCO(1, IBNDS), IREDUNCO(2, IBNDS),
     &                        IBNDS, TOTREDNCO, NRATMS)
            REDUNCO(IBNDS) = DISTAB
C
         ENDIF
C
 20   CONTINUE
C
      DO 30 IANGS = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
 
         IF (IREDUNCO(4, IANGS) .EQ. MONE) THEN
C
C I should check this again, values for 4,ang are -1 or -2(nonlin)
C
            CALL  BULT_LANGCRD(CARTCOORD, BMATRX, IREDUNCO(1, IANGS),
     &                         IREDUNCO(2, IANGS), IREDUNCO(3, IANGS),
     &                         IREDUNCO(4, IANGS), IANGS, NRATMS,
     &                         TOTREDNCO, EPSILON)
cSSS            REDUNCO(IANGS) =
C
         ELSE
C
            CALL  BULT_ANGCRD(CARTCOORD, BMATRX, ANGL,
     &                        IREDUNCO(1, IANGS), IREDUNCO(2, IANGS),
     &                        IREDUNCO(3, IANGS), IANGS, TOTREDNCO,
     &                        NRATMS)
            REDUNCO(IANGS) = ANGL
C
         ENDIF
C
 30   CONTINUE
C
      DO 40 IDIHS = (TOTNOFANG + TOTNOFBND + 1),  TOTREDNCO
C
         CALL  BULT_DIHANGCRD(CARTCOORD, BMATRX, DANG,
     &                        IREDUNCO(1, IDIHS), IREDUNCO(2, IDIHS),
     &                        IREDUNCO(3, IDIHS), IREDUNCO(4, IDIHS),
     &                        IDIHS, TOTREDNCO, NRATMS)
            REDUNCO(IDIHS) = DANG
cSSS            WRITE(6,*) "The Dihedral Angle =", DANG
C    
 40   CONTINUE
C
C Let's write the B-matrix to the JOBARC file. 
C
cSSS     CALL PUTREC(20,'JOBARC', 'BMATRIX ', IINTFP*TOTREDNCO, BMATRX)
C
C Form the G matrix which is required to transfrom the gradient
C from the Cartesian coordinates to redundent internal coordinates.   
C The procedures that are followed here is exactly what is described
C in Pulay et al., 96, 2856, 1992 and Peng et al., 17, 49. 1996.
C Built G = BuB(t) where u is an arbitrary nonsingular matrix, usually
C the unit matrix. The notations are consistent with Pulay et al.
C paper.
C
cSSS      Print*, "The B-Matrix"
cSSS      CALL OUTPUT(BMATRX, 1, TOTREDNCO, 1, 3*NRATMS, TOTREDNCO,
cSSS     &            3*NRATMS, 1)
C
      CALL DCOPY(3*NRATMS*TOTREDNCO, BMATRX, 1, BTMP, 1)
C
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, 3*NRATMS, 1.0D0,
     &           BMATRX, TOTREDNCO, BTMP, TOTREDNCO, 0.0D0, GMATRX,
     &           TOTREDNCO)
C
cSSS      Print*, "The G-Matrix:BB^t"
cSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
cSSS     &            TOTREDNCO, 1)
C---DEBUG
cSSS      CALL ZERO(BMATRX, 3*NRATMS*TOTREDNCO)
cSSS      CALL DCOPY(TOTREDNCO*TOTREDNCO, GMATRX, 1, BMATRX, 1)
C---DEBUG
C
C The intermediate G matrix created above is linear dependent, hence
C there are eigenvalues that are zero. Invert the non-zero digonal
C elements to built the Lambda(-1) matrix.
C
      CALL EIG(GMATRX, EIGVECTORS, 1, TOTREDNCO , 1)
C 
cSSS      Print*, "The Eigenvalues of G-Matrix"
cSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
cSSS     &            TOTREDNCO, 1)
cSSS      Print*, "The Eigenvectors of G-Matrix"
cSSS      CALL OUTPUT(EIGVECTORS, 1, TOTREDNCO, 1, TOTREDNCO, 
cSSS     &            TOTREDNCO, TOTREDNCO, 1)
C
      NULLEVAL = 0
      DO I = 1, TOTREDNCO
         IF (GMATRX(I, I) .LE. EPSILON) THEN
             NULLEVAL = NULLEVL + 1
             GMATRX(I, I) = 0.0D0
         ELSE
             GMATRX(I, I) = 1.0D0/GMATRX(I, I)
         ENDIF
      ENDDO
C  
      NONZERO_TOTREDNCO = TOTREDNCO - NULLEVAL
C
C Built the generalized inverse of G-matrix, what proceeded is generally
C known as singular value decomposition (ineversion of singular matrices)
C G^(-1) = [K L] Lambda^(-1)[K^(t) L^(t)]
C 
      CALL XGEMM('N', 'N', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
     &           EIGVECTORS, TOTREDNCO, GMATRX, TOTREDNCO, 0.0D0,
     &           BTMP, TOTREDNCO)
cSSS      Print*, "The INVERSE OF B-Matrix:1"
cSSS      CALL OUTPUT(BTMP, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
cSSS     &            TOTREDNCO, 1)
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
     &           BTMP, TOTREDNCO, EIGVECTORS, TOTREDNCO, 0.0D0,
     &           GMATRX, TOTREDNCO)
cSSS      Print*, "The INVERSE OF B-Matrix:2"
cSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
cSSS     &            TOTREDNCO, 1)
C---DEBUG
cSSS      CALL ZERO(BTMP, TOTREDNCO*TOTREDNCO)
cSSS      CALL XGEMM('N', 'N', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
cSSS     &           GMATRX, TOTREDNCO, BMATRX, TOTREDNCO, 0.0D0,
cSSS     &           BTMP, TOTREDNCO)
cSSS      CALL XGEMM('N', 'N', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
cSSS     &           BTMP, TOTREDNCO, BMATRX, TOTREDNCO, 0.0D0,
cSSS     &           GMATRX, TOTREDNCO)

cSSS      Print*, "The CHECK INVERSE OF B-Matrix"
cSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
cSSS     &            TOTREDNCO, 1)
C---DEBUG
C
C Built the intermediate matrix needed to built G and A matrices. 
C
      CALL XGEMM('N', 'N', TOTREDNCO, 3*NRATMS, TOTREDNCO, 1.0D0,
     &           GMATRX, TOTREDNCO, BMATRX, TOTREDNCO, 0.0D0,
     &           BTMP, TOTREDNCO)
C  
C Built the the G-MATRIX (G(-1)B)^t; Note the transpose 
C is necessary because the way the transformation is done in 
C CONVQ.f (ACES II). This is the transformation matrix, when
C act on right convert Cartesian gradients to internals, and
C act on left convert internal coordiantes to Cartesian. 
C   
      CALL TRANSP(BTMP, GMATRX, 3*NRATMS, TOTREDNCO)
C
C Make the transpose of the B-Matrix and  pass it out. We 
C need this to transform the Cartesian Hessian matrix to 
C internals.
C
      CALL TRANSP(BMATRX, BTMP, TOTREDNCO, 3*NRATMS)
      CALL DCOPY(3*NRATMS*TOTREDNCO, BTMP, 1, BMATRX, 1)  
C
cSSS      Print*, "The A-Matrix"
cSSS      CALL OUTPUT(GMATRX, 1, 3*NRATMS, 1, TOTREDNCO, 3*NRATMS,
cSSS     &            TOTREDNCO, 1) 
C
      RETURN
      END
