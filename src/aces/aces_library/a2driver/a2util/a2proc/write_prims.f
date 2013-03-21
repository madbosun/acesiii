       Subroutine Write_Prims(NATOMS, NTOT_PRIMS, NBFNS, NAOBFNS,
     &                        COORD, ALPHA, PCOEF, NUNQSHL, NSHL,
     &                        NANGMOMSHL, NCONFUNSHL, NPRIMFUNSHL,
     &                        CNT_COORD, ANGTYPE, IATMCHRG)
C 
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










c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C#include "istart.com"                      

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)

C    
      DIMENSION ALPHA(NTOT_PRIMS), IPOPF(NATOMS), COORD(NATOMS*3), 
     &          NSHL(NATOMS), NANGMOMSHL(NATOMS, NUNQSHL),
     &          NCONFUNSHL(NATOMS, NUNQSHL), 
     &          NPRIMFUNSHL(NATOMS,NUNQSHL), PCOEF(NTOT_PRIMS, NBFNS), 
     &          IANGTYPE(NTOT_PRIMS), ATOM_XYZ(4, MXATMS),
     &          CNT_COORD(3,NTOT_PRIMS),CNTMU(3)
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
C
      DATA XTANG /0.529177249D0/ 
C
       Write(6,*) 
       Write(6,"(a,5I4)"), "Input varaibles:",  Natoms, 
     &         nunqshl, nbfns, NTOT_PRIMS 
       Write(6,*)
       write(*, '(a)'), "The exponents"
       Write(*, '(4(1x,F10.5))') (alpha(i), i=1, NTOT_PRIMS)
       Write(6,*)
       Write(*, '(a)'), "Cartesian Coords."
       Write(*, '(3F10.5)') (Coord(i), i=1, 3*NATOMS)
       Write(6,*)
       Write(*, '(a)'), "The NSHL array"
       Write(*, '(4I5)') (Nshl(i), i=1, NATOMS)
       Write(6,*)
       Write(6,*) "The number of primitive in a shell"
       Write(6, '(4I5)') ((NPRIMFUNSHL(i,j), J=1, Nshl(i)),  
     &                    i=1,NATOMS)
       Write(6,*)
       Write(6,*) "The number of contracted functions in a shell"
       Write(6, '(4I5)') ((NCONFUNSHL(i,j), J=1, Nshl(i)),  
     &                    i=1,NATOMS)
       Write(6,*)
       Write(6,*) "The angular momentum of shells"
       Write(6, '(4I5)') ((NANGMOMSHL(i,j), J=1, Nshl(i)), 
     &                    I=1,NATOMS)
       write(6,*)
       Write(*, '(a)') "The Contractions Coef."
       write(*, '(6F10.5)') ((pcoef(i, j), i=1, NTOT_PRIMS), 
     &                        j=1,nbfns)
      OPEN(UNIT=IUNIT, FILE="BASIS_LOG", FORM="FORMATTED", STATUS="NEW")
C       
      CALL A2GET_ANG_TYPE(NATOMs, NTOT_PRIMS, NUNQSHL, NSHL, NANGMOMSHL, 
     &                    NPRIMFUNSHL, IANGTYPE, COORD, CNT_COORD)
C
      INDEX = 0
      DO LFTPRIM = 1, NTOTPRIM
C
         INDX  = INDX + 1
         ITYPE = IANGTYPE(INDX)
C
         DO IXYZ = 1, 3
            CNTMU(IXYZ) = CNT_COORD(IXYZ, INDX)
         ENDDO
C
         LI = LMNVAL(1, ITYPE)
         MI = LMNVAL(2, ITYPE)
         NI = LMNVAL(3, ITYPE)

      ENDDO 
C

      Return
      END
		
		
	
