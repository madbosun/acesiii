
C THIS ROUTINE DETERMINES THE IRREDUCIBLE REPRESENTATION TO WHICH
C MOLECULAR ORBITALS BELONG, USING A MODEST GENERALIZATION OF THE
C ALGORITHM USED IN JODA FOR DETERMINING SYMMETRIES OF NORMAL MODES.

      SUBROUTINE IRRORB(EVEC,    EVAL,   SCR,    IANGMOM,
     &                  ICENTER, NBAS,           SYOPS,
     &                  SCR1,    IBFATM, ILCATM, Z,
     &                  NATOMS,  iblk,   NBASX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#include "maxbasfn.par"
#include "machsp.com"
#include "mxatms.par"

      CHARACTER*4 IRREP,PTGRP
      CHARACTER*8 IRRSYM(MAXBASFN),SITGRP(MXATMS)
      character*(8*mxatms) szSitGrpTmp

      DIMENSION EVEC(NBASX,NBAS),SCR(1),IANGMOM(NBASX)
      DIMENSION ICENTER(NBASX),SCR1(1),Z(1)
      DIMENSION IBFATM(NATOMS),ILCATM(NATOMS),IPTR(MXATMS*120)
      DIMENSION EVAL(NBAS),SYOPS(1080),ORIENT(9)
      DIMENSION JFLAGS(16)

      COMMON /FLAGS/ IFLAGS(100)

      DATA IONE /1/
      DATA I16 /16/
      DATA DEGTOL /0.00001/

#ifdef _TRACE_VSCF
      write(*,*) '@VSCF: Entering ', __FILE__
#endif /* _TRACE_VSCF */

C PICK UP THE FOLLOWING INFORMATION FROM THE JOBARC FILE:
C
C      1. SCF EIGENVECTORS IN FULL FILLED-OUT AO BASIS.
C      2. ANGULAR MOMENTUM OF EACH BASIS FUNCTION.
C      3. ATOMIC CENTER WHERE EACH BASIS FUNCTION IS LOCATED.
C      4. NUMBER OF ORBITS IN FULL MOLECULAR POINT GROUP.
C      5. ORBIT POPULATION VECTOR IN FULL MOLECULAR POINT GROUP.
C      6. CENSUS VECTOR (TELLS WHICH ATOMS BELONG TO EACH ORBIT).
C      7. PERMUTATION VECTOR FOR ALL SYMMETRY OPERATIONS.
C
C  Let's put in the stuff to do the computational point group as well

      CALL GETREC(20,'JOBARC','JODAFLAG',I16,JFLAGS)
      IGEN=JFLAGS(14)
      do itime = 1, 2

c ----------------------------------------------------------------------

      NBAS2=NBAS*NBASX 

CSSS      if (iblk.eq.1) then
CSSS         CALL GETREC(-1,'JOBARC','EVECAO_A',IINTFP*NBAS2,EVEC)
CSSS      else
CSSS         CALL GETREC(-1,'JOBARC','EVECAO_B',IINTFP*NBAS2,EVEC)
CSSS      end if
CSSS      CALL GETREC(-1,'JOBARC','EVALORDR',IINTFP*NBAS,EVAL)

      CALL GETREC(-1,'JOBARC','ANMOMBF0',NBASX,IANGMOM)
      CALL GETREC(-1,'JOBARC','CNTERBF0',NBASX,ICENTER)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,2I4)") "The NBAS and NBASX: ", Nbas, Nbasx
      Write(6,"(a)")"Eigenvalues, ang. momentum and center of basis fns"
      Write(6,"(6(1x,F10.5))") (Eval(i), i=1, nbas)
      Write(6,*)
      Write(6,"(6I5)") (IANGMOM(i), i=1, NBASX)
      Write(6,*)
      Write(6,"(6I5)") (ICENTER(i), i=1, NBASX)
#endif
      if (itime.eq.1) then
         CALL GETCREC(-1,'JOBARC','FULLPTGP',4,PTGRP)
         CALL GETCREC(-1,'JOBARC','FULLSTGP',8*NATOMS,szSitGrpTmp)
         CALL GETREC(-1,'JOBARC','FULLORDR',IONE,IORDGRP)
         CALL GETREC(-1,'JOBARC','FULLPERM',IORDGRP*NATOMS,IPTR)
         CALL GETREC(-1,'JOBARC','FULLSYOP',IINTFP*9*IORDGRP,SYOPS)
      else
         CALL GETCREC(-1,'JOBARC','COMPPTGP',4,PTGRP)
         CALL GETCREC(-1,'JOBARC','COMPSTGP',8*NATOMS,szSitGrpTmp)
         CALL GETREC(-1,'JOBARC','COMPORDR',IONE,IORDGRP)
         CALL GETREC(-1,'JOBARC','COMPPERM',IORDGRP*NATOMS,IPTR)
         CALL GETREC(-1,'JOBARC','COMPSYOP',IINTFP*9*IORDGRP,SYOPS)
      end if
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)")" Computational,full PG, order, iptr and syops"
      Write(6,"(2a4)") PTGRP, szSitGrpTmp
      Write(6,"(I5)")  IORDGRP
      Write(6,"(6I5)") (IPTR(i), i=1, IORDGRP*NATOMS)
      Write(6,*) 
      Write(6,"(6I5)") (SYOPS(i), i=1, 9*IORDGRP)
#endif
C
      CALL GETREC(-1,'JOBARC','ORIENTMT',IINTFP*9,ORIENT)

c   o 'compress' szSitGrpTmp into the 2-D SITGRP array
      iNdx = 1
      do i = 1, NATOMS
         SITGRP(i)(1:8) = szSitGrpTmp(iNdx:iNdx+7)
         iNdx = iNdx + 8
      end do

      IF  (IFLAGS(1).GE.10) THEN
         write(6,*) '@irrorb: orientmt is '
         call output(orient,1,3,1,3,3,3,1)
         write(6,*) '@irrorb: ops matrices are '
         do i = 1, iordgrp
            call output(syops(1+(i-1)*9),1,3,1,3,3,3,1)
         end do
      END IF

C IF THIS IS A LINEAR MOLECULE, THE FULL POINT GROUP ON DISK
C IS ACTUALLY A SUBGROUP OF THE REAL ONE. RESET THE VALUE HERE.

      IF (PTGRP(1:3).EQ.'C8v'.OR.PTGRP(1:3).EQ.'D8h') THEN
         CALL GETREC(20,'JOBARC','LINEAR  ',IONE,LINEAR)
         IF (LINEAR.EQ.1) PTGRP(2:2)='X'
      END IF

C TRANSFORM ALL SYMMETRY OPERATIONS TO THE COMPUTATIONAL ORIENTATION

      IF (ITIME.EQ.1) CALL TRNOPS(SYOPS,ORIENT,IORDGRP)
      DO I = 1, NBAS
         IRRSYM(I)='XXXX    '
      END DO

C NOW WE HAVE EVERYTHING WE NEED TO DETERMINE ALL OF THE IRREPS.
C LOOP OVER MOLECULAR ORBITALS AND COPY ALL BASIS FUNCTIONS OF
C A PARTICULAR ANGULAR MOMENTUM TYPE INTO A SCRATCH VECTOR

      DO IMO = 1, NBAS
      IF (IRRSYM(IMO)(1:4).EQ.'XXXX') THEN

C CHECK FOR DEGENERACIES SKIP LOOP IF THE IRREP HAS ALREADY BEEN 
C SET. (THIS MEANS THAT IT IS NECESSARILY DEGENERATE WITH THE
C PREVIOUS IRREP).

         NDEG=0
         ITOP=MIN(NBAS,IMO+5)
         DO I=IMO+1,MIN(IMO+5,NBAS+1)
c YAU: WHAT IS THIS NONSENSE? EVAL is only declared with NBAS doubles
c      and yet we are reading EVAL(NBAS+1) when IMO > NBAS-5.
            Z1=ABS(EVAL(I)-EVAL(IMO))
            NDEG=NDEG+1
            IF (Z1.GT.DEGTOL) GOTO 999
         END DO
         NDEG=MAX(NDEG,1)
 999     CONTINUE

         if (itime.eq.2) ndeg=1
         IRREP='XXXX'
         ITOTLN=0
         LENMAX=-1
         IOFFD0=0

         DO IDEG = 1, NDEG
            LANG=-1
            DO IANG = 0, 5
            IF (LANG.EQ.-1) THEN

               DO I = 1, NATOMS
                  ILCATM(I)=-1
                  IBFATM(I)=-1
               END DO
               IOFF=1
               DO IBF = 1, NBASX
                  IF (IANGMOM(IBF).EQ.IANG) THEN
                     SCR(IOFF+IOFFD0)=EVEC(IBF,IMO+IDEG-1)
                     IBFATM(ICENTER(IBF))=IBFATM(ICENTER(IBF))+1
                     IF (ILCATM(ICENTER(IBF)).EQ.-1) THEN
                        ILCATM(ICENTER(IBF))=IOFF
                        IBFATM(ICENTER(IBF))=1
                     END IF
                     IOFF=IOFF+1
                  END IF
               END DO
               LEN=IOFF-1

C CHECK TO MAKE SURE THAT THE EIGENVECTOR ONLY HAS NONZERO
C ENTRIES FOR AN ATOM WHOSE SITE GROUP IS THE SAME AS THE POINT
C GROUP OF THE MOLECULE.  IF THIS IS TRUE, THEN IT IS NECESSARY
C TO USE THE P FUNCTION REPRESENTATION.

               IOK=1
               IF (IANG.EQ.0) THEN
                  IOK=0
                  CALL CHKVEC(SCR(IOFFD0+1),LEN,NATOMS,
     &                        PTGRP,SITGRP,IBFATM,IOK)
               END IF

C CHECK ALSO TO MAKE SURE THAT THIS EIGENVECTOR HAS A NONZERO LENGTH

               X=SNRM2(LEN,SCR(IOFFD0+1),1)
               IF (DABS(X).LT.1.D-06) IOK=0
               IF (IOK.EQ.0) LEN=0
               IF (LEN.NE.0) THEN
                  LANG=IANG
                  LENMAX=MAX(LEN,LENMAX)
                  ITOTLN=ITOTLN+IOFF-1
               END IF

c           END IF (LANG.EQ.-1)
            END IF
c           END DO IANG=0,5
            END DO
            IOFFD0=IOFFD0+LEN
c        END DO IDEG=1,NDEG
         END DO

         IF ((JFLAGS(16).EQ.1.AND.LANG.GE.2).OR.LANG.GE.5) THEN
            IRREP='----'
         ELSE
            CALL IRREPGET(PTGRP,  IORDGRP, Z,      SCR,
     &                    ITOTLN, LEN,     LENMAX, NDEG,
     &                    IPTR,   IBFATM,  ILCATM, SYOPS,
     &                    SCR1,   NATOMS,  IRREP,  NBAS,
     &                    NBASX,  LANG)
         END IF
         IRRSYM(IMO)(1:4)=IRREP
         IF (NDEG.GT.1) THEN
            DO I = 1, NDEG-1
               IRRSYM(IMO+I)(1:4)=IRREP
            END DO
         END IF

c     END IF (IRRSYM(IMO)(1:4).EQ.'XXXX')
      END IF
c     END DO IMO = 1, NBAS
      END DO

      if (iblk.eq.1) then
         if (itime.eq.1) then
           CALL PUTCREC(20,'JOBARC','EVCSYMAF',8*NBAS,IRRSYM)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "The full point group sym. labels: A spin"
         Write(6,"(8A4)") (IRRSYM(i), i=1, Nbas)
         Write(6,*)
#endif  
         else
           CALL PUTCREC(20,'JOBARC','EVCSYMAC',8*NBAS,IRRSYM)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "The comp. point group sym. labels: A spin"
         Write(6,"(8A4)") (IRRSYM(i), i=1, Nbas)
         Write(6,*)
#endif  
         end if
      else
         if (itime.eq.1) then
           CALL PUTCREC(20,'JOBARC','EVCSYMBF',8*NBAS,IRRSYM)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "The full point group sym. labels: B spin" 
         Write(6,"(8A4)") (IRRSYM(i), i=1, Nbas)
         Write(6,*)
#endif  
         else
           CALL PUTCREC(20,'JOBARC','EVCSYMBC',8*NBAS,IRRSYM)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "The comp. point group sym. labels: B spin"
         Write(6,"(8A4)") (IRRSYM(i), i=1, Nbas)
         Write(6,*)
#endif  
         end if
      end if


c ----------------------------------------------------------------------

c     end do itime = 1, 2
      end do

#ifdef _TRACE_VSCF
      write(*,*) '@VSCF: Leaving ', __FILE__
#endif /* _TRACE_VSCF */

      RETURN
      END

