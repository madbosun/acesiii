      subroutine fixrhs(w,d,den,e,ntsiz,pop,vrt,nocc,nvrt,ispin,iuhf)
C
cjp used to prepare from T1 RHS the <j,b|HexpT|0> term for 4.28
cjp NOTE: this routine can be called only with ispin=1,2
cjp
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINCC,CICALC
      INTEGER DIRPRD,POP,VRT
cjp
      integer itrue,iatrue
      DIMENSION W(NTSIZ),E(NOCC+NVRT),
     &          D(NTSIZ),DEN(NTSIZ),POP(8),VRT(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /CORENG/ ELAST
C
      EQUIVALENCE(IFLAGS(2),METHOD)
C
cjp
#include "bwcc.com"
C
C     GET ORBITAL ENERGIES FROM DISK
C
      if(ispin.eq.1.or.iuhf.eq.0)then
       CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*(NOCC + NVRT),
     &             E)
      ELSE
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*(NOCC + NVRT),
     &             E)
      ENDIF
cjp fix diagonal fock elements
      if(isbwcc) then
        do i=1,nocc+nvrt
          e(i)=e(i)+fockcd(i,iref,ispin)
        enddo
      endif
C
C     LOOP OVER IRREPS
C
      IND=0
      INDI=0
      INDA0=0
C
      DO 100 IRREP=1,NIRREP
C
      NOCCI=POP(IRREP)
      NVRTI=VRT(IRREP)
C
cjp
cjp HERE DO NOT identify the internal amplitudes
cjp the following is also completely wrong if symmetry is on
      DO  10 I=1,NOCCI
      INDI=INDI+1
      INDA=INDA0
      DO  10 IA=1,NVRTI
      INDA=INDA+1
      IND=IND+1
           itrue=invhnum(indi,iref,ispin)
           iatrue=invpnum(inda,iref,ispin)
           den(ind)=e(itrue)-e(iatrue)
cjp correction 2: put just LHS of BWCC T1, for internal one keep just the RHS as is
           if(correctiontype.eq.2.and.
     & .not.(isactive(itrue,ispin).and.isactive(iatrue,ispin))) then
           w(ind)=ecorrbw0*d(ind)
           endif
           if(bwgossip.and.correctiontype.eq.2)
     &   write(6,*)'fixrhs:c=2 : ind, ispin, iref,w(ind)'
     &              ,ind, ispin, iref,w(ind)
   10 CONTINUE
      INDA0=INDA0+NVRTI
100   CONTINUE
C
C    compute addition to the rhs of RSCC to get the <...> elements for 4.28
C
cjp use old, internal-zeroed, amplitudes
cjp this is the correct version
#ifndef BAD_ALTERNATIVE
c
       if(correctiontype.le.1) then
cjp
cjp correction type =1 ... build full <> term from current amplitudes and denominators
cjp at convergence scaled to RSCC, this gives zeros except internal excitations
cjp when it is Heff element
cjp
      DO  I=1,NTSIZ
      w(i)=w(i)-d(i)*den(i)
      enddo
c
      endif
cjp correctiontype.eq.2 has been done above already

#else
cjp OBSOLETE
cjp use new, not internal-zeroed, amplitudes
      DO  I=1,NTSIZ
      w(i)=w(i)-w(i)*den(i)/(den(i)+ecorrbw)
        enddo
#endif
C
      RETURN
      END
