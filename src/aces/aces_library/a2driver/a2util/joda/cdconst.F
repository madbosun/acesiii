      SUBROUTINE CDCONST(FREQ,DIDQ,NX)
C
C CALCULATE CENTRIFUGAL DISTORTION TENSOR ELEMENTS
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TMP(3,3),TAU(3,3,3,3),DIDQ(3,3,NX),FREQ(NX)
      DIMENSION CORVEC(3),BSRMBPP(3),BSRMB(3),BSRMBP(3)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/LINEAR/ILINEAR
C
C READ IN MOMENTS OF INERTIA
C
      CALL GETREC(20,'JOBARC','I4CDCALC',IINTFP*9,TMP)
      TMAX=MAX(TMP(1,1),TMP(2,2),TMP(3,3))
      TMIN=MIN(TMP(1,1),TMP(2,2),TMP(3,3))
      DO 5 I=1,3
       IF(TMP(I,I).EQ.TMAX)IC=I
       IF(TMP(I,I).EQ.TMIN)IA=I
5     CONTINUE
C
CWARNING: FANCY LOGIC
C
      IB=6-IC-IA
C   
C
C LOOP OVER INDICES DISREGARDING REDUNDANCIES
C
      FACTOR=60.19968599D0
      FACT=-872720.3477677D0
      FACT2=29979.2458D0
      WRITE(6,1001)
      WRITE(6,1000)
      WRITE(6,1002)
      WRITE(6,1001)
      DO 10 I=1,3-ILINEAR
       DO 20 J=1,3-ILINEAR
        DO 30 K=1,3-ILINEAR
         DO 40 L=1,3-ILINEAR
          Z=0.0D0
          DO 41 IMODE=7-ILINEAR,NX
           Z=Z+DIDQ(I,J,IMODE)*DIDQ(K,L,IMODE)/(FREQ(IMODE)**2)
41        CONTINUE
          if(abs(tmp(i,i)*tmp(j,j)*tmp(k,k)*tmp(l,l)).gt.1.D-8) then
          TAU(I,J,K,L)=fact*Z/(TMP(I,I)*TMP(J,J)*TMP(K,K)*TMP(L,L))
          write(6,'(4i5,e20.10,F20.10)')i,j,k,l,tau(i,j,k,l),
     &             tau(i,j,k,l)*29979.2458D0
          else
           write(6,*) ' skipped ',i,j,k,l
           tau(i,j,k,l)=0.D0
          endif
40       CONTINUE
30      CONTINUE
20     CONTINUE
10    CONTINUE
      WRITE(6,1001)
C
C CALCULATE B''-B FOR ROTATIONAL CONSTANTS.
C
      BXPPMBX=-0.5D0*(TAU(2,2,3,3)+0.5D0*TAU(2,3,2,3)+TAU(1,2,1,2)+
     &      TAU(1,3,1,3))*FACT2
      BYPPMBY=-0.5D0*(TAU(1,1,3,3)+0.5D0*TAU(1,3,1,3)+TAU(2,3,2,3)+
     &      TAU(1,2,1,2))*FACT2
      BZPPMBZ=-0.5D0*(TAU(1,1,2,2)+0.5D0*TAU(1,2,1,2)+TAU(2,3,2,3)+
     &      TAU(1,3,1,3))*FACT2
C
C CALCULATE B''-B' FOR ROTATIONAL CONSTANTS
C
      BXPPMBXP=-0.5D0*(TAU(2,2,3,3)+2.0D0*TAU(2,3,2,3))*FACT2
      BYPPMBYP=-0.5D0*(TAU(1,1,3,3)+2.0D0*TAU(1,3,1,3))*FACT2
      BZPPMBZP=-0.5D0*(TAU(1,1,2,2)+2.0D0*TAU(1,2,1,2))*FACT2
C
C NOTE THAT CODE WHICH FOLLOWS HEREAFTER IS BASED ON THE
C A-B-C REPRESENTATION.
C
C CALCULATE B'-B FOR ROTATIONAL CONSTANTS AS THESE ARE NEEDED FOR
C THE ASYMMETRY PARAMETER
C
      APMA=0.25D0*(3.0D0*TAU(IB,IC,IB,IC)-2.0D0*TAU(IA,IC,IA,IC)-
     &     2.0D0*TAU(IA,IB,IA,IB))
      BPMB=0.25D0*(3.0D0*TAU(IA,IC,IA,IC)-2.0D0*TAU(IB,IC,IB,IC)-
     &     2.0D0*TAU(IA,IB,IA,IB))
      CPMC=0.25D0*(3.0D0*TAU(IA,IB,IA,IB)-2.0D0*TAU(IB,IC,IB,IC)-
     &     2.0D0*TAU(IA,IC,IA,IC))
C
C EVALUATE A',B',C' CONSTANTS FROM RIGID-ROTOR A,B,C (in MHz).
C ALSO CALCULATE SIGMA PARAMETER (WHICH IS AN APPROXIMATION,
C AS IT IS CURRENTLY BASED ON Be' RATHER THAN B0' CONSTANTS)
C
      DO 1010 I=1,3
       IF(ABS(TMP(I,I)).GT.1.D-8) THEN
        TMP(I,I)=FACTOR/TMP(I,I)
       ELSE
        TMP(I,I)=1.D+99
       ENDIF
1010  CONTINUE
      TMP(IA,IA)=TMP(IA,IA)+APMA
      TMP(IB,IB)=TMP(IB,IB)+BPMB
      TMP(IC,IC)=TMP(IC,IC)+CPMC
      IF(ABS(TMP(IB,IB)-TMP(IC,IC)).GT.1.D-08) THEN
       SIGMA=(2.0D0*TMP(IA,IA)-TMP(IB,IB)-TMP(IC,IC))/
     &                       (TMP(IB,IB)-TMP(IC,IC))
      ELSE
       SIGMA=1.D+80
      ENDIF
C
C CALCULATE R5 AND R6 PARAMETERS
C
      R5=-(TAU(IB,IB,IB,IB)-TAU(IC,IC,IC,IC)-2.0D0*TAU(IA,IA,IB,IB)-
     &   4.0D0*TAU(IA,IB,IA,IB)+2.0D0*TAU(IC,IC,IA,IA)+
     &   4.0D0*TAU(IA,IC,IA,IC))/32.0D0
      R6=(TAU(IB,IB,IB,IB)+TAU(IC,IC,IC,IC)-2.0D0*
     &    TAU(IB,IB,IC,IC)-4.0D0*TAU(IB,IC,IB,IC))/64.0D0
c      write(6,*)' R5 and R6 ',r5,r6
C
C CALCULATE BSR-B (S-REDUCED HAMILTONIAN CONSTANTS = BSR)
C
c      ASRMA=(6.0D0*R6-5.0D0*R5/SIGMA)*FACT2
c      BSRMB=(-4.0D0*R6+2.0D0*(2.0D0+1.0D0/SIGMA)*R5)*FACT2
c      CSRMC=(-4.0D0*R6-2.0D0*(2.0D0-1.0D0/SIGMA)*R5)*FACT2
c      write(6,*)asrma,bsrmb,csrmc,'*** bs -> bprime '
      BSRMB(IA)=(6.0D0*R6-5.0D0*R5/SIGMA+0.25D0*
     &      (3.0D0*TAU(IB,IC,IB,IC)-2.0D0*TAU(IA,IC,IA,IC)
     &      -2.0D0*TAU(IA,IB,IA,IB)))*FACT2
      BSRMB(IB)=(-4.0D0*R6+2.0D0*(2.0D0+1.0D0/SIGMA)*R5+0.25D0*
     &      (3.0D0*TAU(IA,IC,IA,IC)-2.0D0*TAU(IB,IC,IB,IC)
     &      -2.0D0*TAU(IA,IB,IA,IB)))*FACT2
      BSRMB(IC)=(-4.0D0*R6-2.0D0*(2.0D0-1.0D0/SIGMA)*R5+0.25D0*
     &      (3.0D0*TAU(IA,IB,IA,IB)-2.0D0*TAU(IB,IC,IB,IC)
     &      -2.0D0*TAU(IA,IC,IA,IC)))*FACT2
c      write(6,*)bsrmb,'*** bs -> b '
C
C CALCULATE FANCY PARAMETERS FOR WATSON'S A-REDUCED HAMILTONIAN
C
      DDJ=-(3.0D0*(TAU(IC,IC,IC,IC)+TAU(IB,IB,IB,IB))+2.0D0*
     &     (TAU(IB,IB,IC,IC)+2.0D0*TAU(IB,IC,IB,IC)))/32.0D0
      DDK =DDJ-0.25D0*(TAU(IA,IA,IA,IA)-TAU(IA,IA,IB,IB)-2.0D0*
     &        TAU(IA,IB,IA,IB)-TAU(IA,IA,IC,IC)-2.0D0*
     &        TAU(IA,IC,IA,IC))
      DDJK=-DDJ-DDK-0.25D0*TAU(IA,IA,IA,IA)
c      write(6,*)ddj,ddk,tau(ia,ia,ia,ia)
      DELJ=DDJ-2.D0*R6
      DELK=DDK-10.0D0*R6
      DELJK=DDJK+12.0D0*R6
      DEL2J=-(TAU(IB,IB,IB,IB)-TAU(IC,IC,IC,IC))/16.0D0
      DEL2K=-2.0D0*R5-4.0D0*R6*SIGMA
c      write(6,*)' A-reduced Hamiltonian centrifugal constants'
c      write(6,*)' r6 is ',r6,2.0d0*r5,4.0d0*r6*sigma
c      write(6,*)' DJ is ',ddj,ddj*fact2
c      write(6,*)' DK is ',ddk,ddk*fact2
c      write(6,*)' DJK is ',ddjk,ddjk*fact2
c      write(6,*)' DELJ is ',delj,delj*fact2
c      write(6,*)' DELK is ',delk,delk*fact2
c      write(6,*)' DELJK is ',d(6,*)' DELK is ',delk,delk*fact2
c      write(6,*)' DELJK is ',deljk,deljk*fact2
c      write(6,*)' delJ is ',del2j,del2j*fact2
c      write(6,*)' delK is ',del2k,del2k*fact2
c      write(6,*)' SIGMA is ',sigma
C
C WRITE OUT THE STUFF
C
      WRITE(6,304)
      WRITE(6,305)
      WRITE(6,304)
      WRITE(6,*)
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'R6',r6*fact2,r6
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'R5',R5*FACT2,R5
      WRITE(6,'(T12,A2,T42,E12.6)')'SIGMA',SIGMA
      WRITE(6,*)
      WRITE(6,*)' A-reduced centrifugal distortion parameters '
      WRITE(6,*)
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'DJ',DDJ*FACT2,DDJ
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'DK',DDK*FACT2,DDK
      WRITE(6,'(T12,A3,T30,E12.6,T53,E12.6)')'DJK',DDJK*FACT2,DDJK
      WRITE(6,'(T12,A4,T30,E12.6,T53,E12.6)')'DELJ',DELJ*FACT2,DELJ
      WRITE(6,'(T12,A4,T30,E12.6,T53,E12.6)')'DELK',DELK*FACT2,DELK
      WRITE(6,'(T12,A5,T30,E12.6,T53,E12.6)')'DELJK',DELJK*FACT2,DELJK
      WRITE(6,'(T12,A4,T30,E12.6,T53,E12.6)')'delJ',DEL2J*FACT2,DEL2J
      WRITE(6,'(T12,A4,T30,E12.6,T53,E12.6)')'delK',DEL2K*FACT2,DEL2K
      WRITE(6,*)
      WRITE(6,*)' S-reduced centrifugal distortion parameters '
      WRITE(6,*)
C
C CALCULATE FANCY PARAMETERS FOR WATSON'S S-REDUCED HAMILTONIAN
C
      DDJ=DDJ+R5/SIGMA
      DDK=DDK+R5*5.0D0/SIGMA
      DDJK=DDJK-6.0D0*R5/SIGMA
      D1=-DEL2J
      D2=R6+R5/(2.0D0*SIGMA)
C
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'DJ',DDJ*FACT2,DDJ
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'DK',DDK*FACT2,DDK
      WRITE(6,'(T12,A3,T30,E12.6,T53,E12.6)')'DJK',DDJK*FACT2,DDJK      
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'D1',D1*FACT2,D1
      WRITE(6,'(T12,A2,T30,E12.6,T53,E12.6)')'D2',D2*FACT2,D2
      WRITE(6,304)

c      write(6,*)' d2 is ',d2
c      write(6,*)' S-reduced Hamiltonian centrifugal constants'
c      write(6,*)' r6 is ',r6,r6*fact2
c      write(6,*)' DJ is ',ddj,ddj*fact2
c      write(6,*)' DK is ',ddk,ddk*fact2
c      write(6,*)' DJK is ',ddjk,ddjk*fact2
c      write(6,*)' d1 is ',D1,D1*fact2
c      write(6,*)' d2 is ',D2,D2*fact2
c      write(6,*)' SIGMA is ',sigma
c
c      write(6,*)' debug stuff '
      BSRMBPP(IB)=-(2.0D0*ddj+ddjk+2.0d0*(d1+2.0d0*d2))*FACT2
      BSRMBPP(IC)=-(2.0D0*ddj+ddjk+2.0d0*(-d1+2.0d0*d2))*FACT2
      BSRMBPP(IA)=-(2.0D0*ddj+6.0d0*d2)*FACT2
c      write(6,*)z3*fact2,z1*fact2,z2*fact2,' **** bs -> bdprime '
      BSRMBP(IB)=-(4.0D0*r6-4.0D0*r5-2.0d0*r5/sigma)*FACT2
      BSRMBP(IC)=-(4.0D0*r6+4.0d0*r5-2.d0*r5/sigma)*FACT2
      BSRMBP(IA)=-(6.0d0*r6-5.0d0*r5/sigma)*FACT2
C
      WRITE(6,*)' Centrifugal corrections for various types of ',
     &          'rotational constants (MHz) ...'
      WRITE(6,300)BXPPMBX,BYPPMBY,BZPPMBZ
300   FORMAT(T3,' Bx''''-Bx ',F20.10,/,T3,' By''''-By ',F20.10,/,
     &       T3,' Bz''''-Bz ',F20.10)
      WRITE(6,301)BXPPMBXP,BYPPMBYP,BZPPMBZP
301   FORMAT(T3,' Bx''''-Bx'' ',F20.10,/,T3,' By''''-By'' ',F20.10,/,
     &       T3,' Bz''''-Bz'' ',F20.10)
      WRITE(6,302)BXPPMBX-BXPPMBXP,BYPPMBY-BYPPMBYP,BZPPMBZ-BZPPMBZP
302   FORMAT(T3,' Bx''-Bx ',F20.10,/,T3,' By''-By ',F20.10,/,
     &       T3,' Bz''-Bz ',F20.10)
      WRITE(6,303)BSRMB(1),BSRMB(2),BSRMB(3)
303   FORMAT(T3,' BxSR-Bx ',F20.10,/,T3,' BySR-By ',F20.10,/,
     &       T3,' BzSR-Bz ',F20.10)
      WRITE(6,308)BSRMBPP(1),BSRMBPP(2),BSRMBPP(3)
308   FORMAT(T3,' BxSR-Bx'''' ',F20.10,/,T3,' BySR-By'''' ',F20.10,/,
     &       T3,' BzSR-Bz'''' ',F20.10)
      WRITE(6,309)BSRMBP(1),BSRMBP(2),BSRMBP(3)
309   FORMAT(T3,' BxSR-Bx'' ',F20.10,/,T3,' BySR-By'' ',F20.10,/,
     &       T3,' BzSR-Bz'' ',F20.10)
C
305   FORMAT(T9,'Parameter',T34,'(MHz)',T56,'(CM-1)')
304   FORMAT(60('-'))
C
1000  FORMAT(T10,' Quartic centrifugal distortion parameters ')
1002  FORMAT(T33,'CM-1',T53,'MHz')
1001  FORMAT(75('-'))
      RETURN
      END
