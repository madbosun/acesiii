













































































































































































































































































































































































































































































































































      SUBROUTINE A2MAKE_WFN(DENSA, DENSB, NATORBA, NATORBB, EVECA,
     &                      EVECB, EIGVALA, EIGVALB, SCR1, SCR2,
     &                      OCCA, OCCB, NOCCA, NOCCB, MAXOCCA,
     &                      MAXOCCB, NAOBFNS, NBAS, IUHF, TENERGY)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
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



c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end







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




C
      DOUBLE PRECISION NATORBA, NATORBB, OCCA(NBAS), OCCB(NBAS)
      LOGICAL SCF, COR
      DIMENSION DENSA(NBAS*NBAS), DENSB(NBAS*NBAS), NOCC(16),
     &          NATORBA(NAOBFNS, NAOBFNS), NATORBB(NAOBFNS, NAOBFNS),
     &          EVECA(NBAS, NBAS), EVECB(NBAS, NBAS),
     &          EIGVALA(NBAS), EIGVALB(NBAS), 
     &          SCR1(NBAS*NBAS), SCR2(NBAS*NBAS)
C
      DATA ZILCH, ONE, HALF/ 0.0D0, 1.0D0, 0.50D0/
C
      SCF = (iflags(2) .eq.0)
      COR = (iflags(2) .gt.0)
C
      
      IF (SCF) THEN

         CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP, EVECA)
         CALL GETREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP, EVALA)
         DO IOCC = 1, NOCCA
            OCCA(IOCC) = 2.0D0
         ENDDO
         MAXOCCA = NOCCA
C
         CALL GETREC(20, "JOBARC", "KINETINT", NBAS*NBAS*IINTFP, SCR1)
         CALL GETREC(20, "JOBARC", "SCFDENSA", NBAS*NBAS*IINTFP, SCR2)
         TENERGYA = DDOT(NBAS*NBAS, SCR1, 1, SCR2, 1)
C
         IF (IUHF .EQ. 1) THEN
            CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,
     &                  EVECB)
            CALL GETREC(20,'JOBARC','SCFEVLB0',NBAS*IINTFP, EVALB)
            Do IOCC = 1, NOCCA
               OCCA(IOCC) = 1.0D0
            Enddo
            Do IOCC = 1, NOCCB
                OCCB(IOCC) = 1.0D0
            ENDDO
            MAXOCCB = NOCCB
C
         CALL GETREC(20, "JOBARC", "SCFDENSB", NBAS*NBAS*IINTFP, SCR2)
         TENERGYB = DDOT(NBAS*NBAS, SCR1, 1, SCR2, 1)
C
         ENDIF
         TENERGY = TENERGYA + TENERGYB
      Write(6,*)
      Write(6,"(a,F10.5)") "The total kinetic energy = ", TENERGY
C
      ELSE IF (COR) THEN
C
         CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,DENSA)
         CALL DCOPY(NBAS*NBAS, DENSA, 1, SCR2, 1)
         CALL MO2AO3(SCR2, DENSA, EVECA, SCR1, NBAS, NBAS, 1, 2)
         CALL GETREC(20, "JOBARC", "KINETINT", NBAS*NBAS*IINTFP, SCR1)
         TENERGYA = DDOT(NBAS*NBAS, DENSA, 1, SCR1, 1)
C
         CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,DENSA)
C
C
         CALL EIG(DENSA, SCR1, 1, NBAS, -1)
C
         CALL DCOPY(NBAS, DENSA, NBAS+1, OCCA, 1)
         CALL DZERO(EVALA, NBAS)
         MAXOCCA = NBAS
         CALL GETREC(20,'JOBARC','SCFEVECA',NBAS*NBAS*IINTFP,DENSA)
         CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,DENSA,NBAS,SCR1, 
     &               NBAS,ZILCH,EVECA,NBAS)
C
      Write(6,*) "The Alpha Natural orbitals AOxMO"
      call output(EVECA, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
C
        IF (IUHF.EQ.1) THEN
           CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,DENSB)
           CALL DCOPY(NBAS*NBAS, DENSB, 1, SCR2, 1)
           CALL MO2AO3(SCR2, DENSB, EVECA, SCR1, NBAS, NBAS, 2, 2)
           CALL GETREC(20, "JOBARC", "KINETINT", NBAS*NBAS*IINTFP, SCR1)
           TENERGYB = DDOT(NBAS*NBAS, DENSB, 1, SCR1, 1)
C
           CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,DENSB)
C
       Write(6,*) "The Beta Density read"
       call output(DENSB, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
           CALL EIG(DENSB, SCR1, 1, NBAS, -1)
           CALL DCOPY(NBAS, DENSB, NBAS+1, OCCB, 1)
           CALL DZERO(EVALB, NBAS)
           MAXOCCB = NBAS
           CALL GETREC(20,'JOBARC','SCFEVECB',NBAS*NBAS*IINTFP,DENSB)
           CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,DENSB,NBAS,SCR1,
     &                 NBAS,ZILCH,EVECB,NBAS)
      Write(6,*) "The Beta Natural orbitals AOxMO"
      call output(EVECB, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
        ENDIF
        TENERGY = TENERGYA + TENERGYB

      Write(6,*)
      Write(6,"(a,F10.5)") "The total kinetic energy = ", TENERGY
C
      ENDIF
C
      CALL GETREC(20,'JOBARC','CMP2ZMAT', NBAS*NAOBFNS*IINTFP, SCR1)
      CALl XGEMM('N','N', NAOBFNS, NBAS, NBAS, ONE, SCR1, NAOBFNS, 
     &            EVECA, NBAS, ZILCH, NATORBA, NAOBFNS)
C
      Write(6,*)
      Write(6,*) "The Alpha Natural orbitals NAOBFNSxNMO"
      call output(NATORBA, 1, NAOBFNS, 1, NBAS, NAOBFNS,  NBAS, 1)
C
      IF (IUHF .EQ. 1) THEN
         CALL XGEMM('N','N', NAOBFNS, NBAS, NBAS, ONE, SCR1, NAOBFNS,
     &               EVECB, NBAS, ZILCH, NATORBB, NAOBFNS)
      ENDIF
  
      RETURN
      END
