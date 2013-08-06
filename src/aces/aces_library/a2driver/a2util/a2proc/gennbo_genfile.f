
      SUBROUTINE GENNBO_GENFILE(NATOMS,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end


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











      DIMENSION COORD(3,Mxatms),NUCCHG(Mxatms),OVLAP(MAXBASFN*MAXBASFN)
      CHARACTER*80 FNAME
      CHARACTER*4 ATMLBL(MAXBASFN),LABEL(10),ANGLBL(MAXBASFN)
      INTEGER CENTERBF(MAXBASFN)

      DATA LABEL/'S','X','Y','Z','XX','YY','ZZ','XY','XZ','YZ'/

      OPEN(UNIT=47,FILE='FILE47')

      WRITE(47,*)'$GENNBO'
      WRITE(47,'(T2,A,I3)')'NATOMS=',NATOMS
      WRITE(47,'(T2,A,I3)')'NBAS=',NBAS
      WRITE(47,*)'BODM'
      WRITE(47,*)'BOHR'
      WRITE(47,*)'$END'
      WRITE(47,*)
      WRITE(47,*)'$NBO'
      WRITE(47,*)'PRINT=4'
      WRITE(47,*)'NLMO'
      WRITE(47,*)'$END'
      WRITE(47,*)


      WRITE(47,*)'$COORD'
      WRITE(47,*)'NATURAL ORBITAL ANALYSIS USING ACES2 INFORMATION'
      CALL GETREC(-1,'JOBARC','COORD',3*NATOMS,COORD)
      CALL GETREC(-1,'JOBARC','ATOMCHRG',NATOMS,NUCCHG)

      DO I = 1, NATOMS
         WRITE(47,*)NUCCHG(I),NUCCHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
      END DO
      WRITE(47,*)'$END'


      CALL GETREC(-1,'JOBARC','CENTERBF',NBAS,CENTERBF)
      WRITE(47,*)
      WRITE(47,*)'$BASIS'
      WRITE(47,*)'CENTER ='

      DO I = 1, NBAS
         WRITE(47,*)CENTERBF(I)
      END DO

      CALL GFNAME('IIII',FNAME,ILENGTH)
      OPEN(UNIT=10,FILE='IIII',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      CALL LOCATE(10,'LABBASIS')

      DO I = 1, NBAS
         READ(10)J,ATMLBL(I),ANGLBL(I),IJUNK
      END DO

      WRITE(47,*)
      WRITE(47,*)'LABEL ='
      IDTYPE=0

      DO I = 1, NBAS
         IF(ANGLBL(I).EQ.LABEL(1))THEN
            WRITE(47,*)'1'
         ELSE IF(ANGLBL(I).EQ.LABEL(2))THEN
            WRITE(47,*)'101'
         ELSE IF(ANGLBL(I).EQ.LABEL(3))THEN
            WRITE(47,*)'102'
         ELSE IF(ANGLBL(I).EQ.LABEL(4))THEN
            WRITE(47,*)'103'
         ELSE IF(ANGLBL(I).EQ.LABEL(5))THEN
            WRITE(47,*)'201'
            IDTYPE=IDTYPE+1
         ELSE IF(ANGLBL(I).EQ.LABEL(6))THEN
            WRITE(47,*)'204'
            IDTYPE=IDTYPE+1
         ELSE IF(ANGLBL(I).EQ.LABEL(7))THEN
            WRITE(47,*)'206'
            IDTYPE=IDTYPE+1
         ELSE IF(ANGLBL(I).EQ.LABEL(8))THEN
            WRITE(47,*)'202'
         ELSE IF(ANGLBL(I).EQ.LABEL(9))THEN
            WRITE(47,*)'203'
         ELSE IF(ANGLBL(I).EQ.LABEL(10))THEN
            WRITE(47,*)'205'
         END IF
      END DO
      WRITE(47,*)'$END'
      WRITE(*,*)'there are',IDTYPE,'D ORBITALS'
      WRITE(47,*)


c     The contract datalist is not necessary since it is used only for plotting data.



      CALL GETREC(-1,'JOBARC','AOOVRLAP',NBAS*NBAS,OVLAP)
      WRITE(47,*)'$OVERLAP'
      WRITE(47,*)OVLAP
      WRITE(47,*)'$END'


      CALL GETREC(-1,'JOBARC','SCFDENSA',NBAS*NBAS,OVLAP)
      WRITE(47,*)
      WRITE(47,*)'$DENSITY'
      WRITE(47,*)OVLAP
      WRITE(47,*)'$END'

      CALL GETREC(-1,'JOBARC','FOCKA',NBAS*NBAS,OVLAP)
      WRITE(47,*)
      WRITE(47,*)'$FOCK'
      WRITE(47,*)OVLAP
      WRITE(47,*)'$END'

      CALL GETREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS,OVLAP)
      WRITE(47,*)
      WRITE(47,*)'$LCAOMO'
      WRITE(47,*)OVLAP
      WRITE(47,*)'$END'


      CLOSE(47)

      RETURN
      END

