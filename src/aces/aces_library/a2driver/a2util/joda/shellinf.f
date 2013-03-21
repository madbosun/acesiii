      SUBROUTINE SHELLINF(NORBIT,NTOTATOM,NTOTSHEL,INUMSHEL,ISHLOFF,
     &                   ISHELTP,ISHELSZ)
C
C GENERATE SHELL INFORMATION FOR ATOMS AS ORDERED IN THE ZMAT
C FILE.
C
CJDW 10/28/96. Modifications to check value of IMOL and to check we do
C              not exceed fixed dimensions. For dummy atoms IMOL is set
C              to 999 and we have troubles. Thanks to Roger Edberg of ANU
C              for finding this error.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION INUMSHEL(NORBIT),ISHLOFF(NORBIT)
      DIMENSION ISHELTP(NTOTSHEL),ISHELSZ(NTOTSHEL)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

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
      DIMENSION IZMT2MOL(mxatms)
C
C I am going to set these two arrays to the maximum number of
C atoms times maximum number of shells. Ajith Perera, 11/07.
C
      DIMENSION ISHELTP2(mxatms*mxshel),ISHELSZ2(mxatms*mxshel)
C
      IF(NTOTATOM.GT.200)THEN
       WRITE(6,1000)
       CALL ERREX
      ENDIF
C
      CALL GETREC(20,'JOBARC','ZMAT2MOL',NTOTATOM,IZMT2MOL)
      Write(6,*) norbit
      Print*, "ZMAT2MOL", (IZMT2MOL(I), I=1, NTOTATOM)
      Print*, "ISHELTP", (ISHELTP(I), I=1, NTOTSHEL)
      Print*, "ISHELSZ", (ISHELSZ(I), I=1, NTOTSHEL)
C
C LOOP OVER ATOMS IN ZMAT
C
      IOFFZMAT=1
      DO 10 IZMAT=1,NTOTATOM
C
C GET MOL FILE POSITION
C
       IMOL=IZMT2MOL(IZMAT)
C
       IF(IMOL.GE.1 .AND. IMOL.LE.NORBIT)THEN
C
C GET SHELL INFORMATION FOR THIS ATOM
C
        NSHELL=INUMSHEL(IMOL)
        IOFFMOL=ISHLOFF(IMOL)
        if (IOFFZMAT+NSHELL.gt.500) then
           print *, '@SHELLINF: Assertion failed.'
           print *, '           maximum number of shells = 500'
           print *, '           require at least ',IOFFZMAT+NSHELL
           call errex
        end if
        CALL ICOPY(NSHELL,ISHELTP(IOFFMOL),1,ISHELTP2(IOFFZMAT),1)
        CALL ICOPY(NSHELL,ISHELSZ(IOFFMOL),1,ISHELSZ2(IOFFZMAT),1)
        IOFFZMAT=IOFFZMAT+NSHELL
C
       ELSE
C
        IF(IMOL.NE.999)THEN
         WRITE(6,1010) IMOL
         CALL ERREX
        ENDIF
C
       ENDIF
C
10    CONTINUE
C
      NSHELTOT=IOFFZMAT-1
C
      IF(NSHELTOT.GT.500)THEN
       WRITE(6,1020) NSHELTOT
       CALL ERREX
      ENDIF
C
C
      CALL PUTREC(20,'JOBARC','FULSHLNM',1,NSHELTOT)
      CALL PUTREC(20,'JOBARC','FULSHLTP',NSHELTOT,ISHELTP2)
      CALL PUTREC(20,'JOBARC','FULSHLSZ',NSHELTOT,ISHELSZ2)
C
C
      RETURN
 1000 FORMAT(' @SHELLINF-F, Too many atoms (over 200) ',I10)
 1010 FORMAT(' @SHELLINF-F, Invalid value of IMOL ',I10)
 1020 FORMAT(' @SHELLINF-F, Too many shells (over 500) ',I10)
      END
