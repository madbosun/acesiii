      Subroutine projec_IFC(HESS, PHESS, PMAT, GMATRX_N, GMATRX_M,
     &                      GRD, NXM6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end









C
      Dimension HESS(NXM6, NXM6), PMAT(NXM6, NXM6), PHESS(NXM6, NXM6),
     &          GMATRX_M(NXM6, NXM6), GMATRX_N(NXM6, NXM6),
     &          GRD(NXM6)
 
      LENGMAT=NXM6*NXM6
      CALL GETREC(20,'JOBARC','GI-MATRX',LENGMAT*IINTFP,GMATRX_M)
      CALL GETREC(20,'JOBARC','G-MATRX ',LENGMAT*IINTFP,GMATRX_N)

      Call XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, GMATRX_M, 
     &            NXM6, GMATRX_N, NXM6, 0.0D0, PMAT, NXM6)
C
C
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)

      RETURN
      END
     
