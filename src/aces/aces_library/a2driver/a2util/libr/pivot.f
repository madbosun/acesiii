      SUBROUTINE PIVOT(AB,N,ND,OUTER,SCRATCH,EPS)
C
C     This subroutine switches two columns of a matrix to get
C         a nonzero entry in the diagonal.
C     Martin J. McBride.  12/04/85.
C     General Electric CRD, Information System Operation.
C
CEND
      INTEGER N,ND,COL,OUTER,I
      DOUBLE PRECISION AB(ND,1),SCRATCH(1),TEMP,EPS

C  Get first column with non-zero element in row OUTER.
      COL = OUTER + 1
   10 IF (COL .GT. N) GO TO 90
      IF (ABS(AB(OUTER,COL)) .GT. EPS) GO TO 20
         COL = COL + 1
         GO TO 10

C  Switch column OUTER with column COL, which has non-zero element in
C  row OUTER.
   20 DO 30 I = 1,N
         TEMP = AB(I,OUTER)
         AB(I,OUTER) = AB(I,COL)
         AB(I,COL) = TEMP
   30 CONTINUE
      TEMP = SCRATCH(N+OUTER)
      SCRATCH(N+OUTER) = SCRATCH(N+COL)
      SCRATCH(N+COL) = TEMP

   90 CONTINUE
      RETURN
      END
