
      SUBROUTINE EXTRAP_MAIN(ARGS,DIMARGS)
      IMPLICIT NONE

      INTEGER DIMARGS
      CHARACTER*80 ARGS(DIMARGS)

      INTEGER NATOMS, zFREE, zSCFGD, zTOTGD



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






c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end

      IF (DIMARGS.EQ.0) THEN
         PRINT *, '@EXTRAP_MAIN: unspecified argument'
         CALL ERREX
      END IF

      CALL GETREC(1, "JOBARC", "NREALATM", 1, NATOMS)

      zSCFGD = I0
      zTOTGD = zSCFGD + 3*NATOMS*IINTFP
      zFREE  = zTOTGD + 3*NATOMS*IINTFP
      IF (zFREE-I0.GT.ICRSIZ) cALL INSMEM('EXTRAP_MAIN',zFREE-1,ICRSIZ)

      CALL EXTRAP_GRAD(ICORE(zSCFGD),ICORE(zTOTGD),NATOMS,ARGS(1))

      RETURN
      END

