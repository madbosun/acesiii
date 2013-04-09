      SUBROUTINE BUILT_BGMTRX(CARTCOORD, REDUNCO, IREDUNCO,
     &                        TOTREDNCO, TOTNOFBND, TOTNOFANG,
     &                        TOTNOFDIH, NRATMS, 
     &                        BMATRX, GMATRX, EPSILON)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=100, MAXCNTVS = 10, MAXREDUNCO = 100)



c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

c io_units.par : end
/* _IO_UNITS_PAR_ */

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
      WRITE(6,*) TOTREDNCO, NRATMS
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
     &                         IREDUNCO(2, IANGS), IANGS, NRATMS,
     &                         TOTREDNCO, EPSILON)
C$$$            REDUNCO(IANGS) =
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
            WRITE(6,*) "The Dihedral Angle =", DANG
C    
 40   CONTINUE
C
C Let's write the B-matrix to the JOBARC file. 
C
C$$$     CALL PUTREC(20,'JOBARC', 'BMATRIX ', IINTFP*TOTREDNCO, BMATRX)
C
C Form the G matrix which is required to transfrom the gradient
C from the Cartesian coordinates to redundent internal coordinates.   
C The procedures that are followed here is exactly what is described
C in Pulay et al., 96, 2856, 1992 and Peng et al., 17, 49. 1996.
C Built G = BuB(t) where u is an arbitrary nonsingular matrix, usually
C the unit matrix. The notations are consistent with Pulay et al.
C paper.
C
C$$$      Print*, "The B-Matrix"
C$$$      CALL OUTPUT(BMATRX, 1, TOTREDNCO, 1, 3*NRATMS, TOTREDNCO,
C$$$     &            3*NRATMS, 1)
C
      CALL DCOPY(3*NRATMS*TOTREDNCO, BMATRX, 1, BTMP, 1)
C
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, 3*NRATMS, 1.0D0,
     &           BMATRX, TOTREDNCO, BTMP, TOTREDNCO, 0.0D0, GMATRX,
     &           TOTREDNCO)
C
      Print*, "The G-Matrix:BB^t"
      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
     &            TOTREDNCO, 1)
C---DEBUG
      CALL DCOPY(3*NRATMS*TOTREDNCO, GMATRX, 1, BMATRX, 1) 
C---DEBUG
C
C The intermediate G matrix created above is linear dependent, hence
C there are eigenvalues that are zero. Invert the non-zero digonal
C elements to built the Lambda(-1) matrix.
C
      CALL EIG(GMATRX, EIGVECTORS, 1, TOTREDNCO , 1)
C 
      Print*, "The Eigenvalues of G-Matrix"
      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
     &            TOTREDNCO, 1)
      Print*, "The Eigenvectors of G-Matrix"
      CALL OUTPUT(EIGVECTORS, 1, TOTREDNCO, 1, TOTREDNCO, 
     &            TOTREDNCO, TOTREDNCO, 1)
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
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
     &           BTMP, TOTREDNCO, EIGVECTORS, TOTREDNCO, 0.0D0,
     &           GMATRX, TOTREDNCO)
C---DEBUG
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
     &           BMATRX, TOTREDNCO, GMATRX, TOTREDNCO, 0.0D0,
     &           BTMP, TOTREDNCO)
      Print*, "The INVERSE OF B-Matrix"
      CALL OUTPUT(BTMP, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
     &            TOTREDNCO, 1)
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
C$$$      Print*, "The A-Matrix"
C$$$      CALL OUTPUT(GMATRX, 1, 3*NRATMS, 1, TOTREDNCO, 3*NRATMS,
C$$$     &            TOTREDNCO, 1) 
C
      RETURN
      END
