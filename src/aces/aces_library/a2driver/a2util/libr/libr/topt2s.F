      SUBROUTINE TOPT2S(ISPIN,ILIST,NTOP,NT2SIZ,NOCCA,NOCCB,
     &                  NVRTA,NVRTB,NSMSZ1,NSMSZ2,T,TOPT2,
     &                  SYVEC1,SYVEC2,ITOPT2,I,J,A,B,TYPE)
C
C THIS ROUTINE PICKS UP A T2 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT1 AND THEIR ASSOCIATED INDICES
C  IN THE I AND A VECTORS.
C
C INPUT:
C       ISPIN - THE LIST SUBTYPE
C       ILIST - THE LIST NUMBER
C       NTOP  - THE NUMBER OF LARGEST AMPLITUDES (BY ABSOLUTE VALUE)
C               WHICH ARE TO BE FOUND
C       NT2SIZ- THE TOTAL SIZE OF THE T2 VECTOR
C       NOCCA - THE NUMBER OF ALPHA OCCUPIED ORBITALS
C       NOCCB - THE NUMBER OF BETA OCCUPIED ORBITALS
C       NVRTA - THE NUMBER OF ALPHA VIRTUAL ORBITALS
C       NVRTB - THE NUMBER OF BETA VIRTUAL ORBITALS
C       NSMSZ1- THE TOTAL SIZE OF THE A,B SYMMETRY VECTOR
C               FOR THIS SPIN CASE (NVRTA*(NVRTA-1))/2 FOR
C               ISPIN=1, NVRTA*NVRTB FOR ISPIN=3, ETC.)
C       NSMSZ2- THE TOTAL SIZE OF THE I,J SYMMETRY VECTOR
C               FOR THIS SPIN CASE (NOCCA*(NOCCA-1))/2 FOR
C               ISPIN=1, NOCCA*NOCCB FOR ISPIN=3, ETC.)
C       SYVEC1- THE A,B SYMMETRY VECTOR
C       SYVEC2- THE I,J SYMMETRY VECTOR
C       TYPE  - ???????
C
C OUTPUT:
C       TOPT2 - THE NTOP LARGEST ELEMENTS IN THE T2 VECTOR
C               SORTED BY ABSOLUTE VALUE
C       I     - THE I INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       J     - THE J INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       A     - THE A INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       B     - THE B INDICES CORRESPONDING TO THE VALUES IN TOPT1
C
C SCRATCH:
C
C       T     - USED TO HOLD THE SYMMETRY-PACKED T VECTOR
C       ITOPT2- USED TO HOLD THE OFFSETS CORRESPONDING TO THE
C               ELEMENTS IN TOPT2
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A(NTOP),B(NTOP),SYVEC1(NSMSZ1),SYVEC2(NSMSZ2)
cjp
      integer itrue, jtrue,atrue,btrue
      CHARACTER*2 SPCASE(3)
      CHARACTER*1 TYPE
      DIMENSION T(NT2SIZ),TOPT2(NTOP),ITOPT2(NTOP),I(NTOP),
     &          J(NTOP)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA SPCASE /'AA','BB','AB'/
cjp
#include "bwcc.com"
C
C THIS ROUTINE PICKS UP A T2 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT2 AND THEIR ASSOCIATED INDICES
C  IN THE I,J,A AND B VECTORS.
C
C
C STATEMENT FUNCTIONS
C
      ILRG(INTX)=INT(0.5D0*(1.D0+SQRT(8.0*INTX-7)))+1
      NNM1O2(IX)=(IX*(IX-1))/2
C
      IF(NT2SIZ.EQ.0)RETURN
C
      CALL GETALL(T,NT2SIZ,1,ILIST)
C
C COMPUTE NORM OF T2 VECTOR
C
      X=SNRM2(NT2SIZ,T,1)
      T2NORM=X
C
C PUT NTOP LARGEST ELEMENTS INTO TOPT2 AND THEIR ASSOCIATED INDICES INTO
C   ITOPT2.
C
      CALL SCANVC(T,TOPT2,ITOPT2,NTOP,NT2SIZ)
C
C BRANCH HERE FOR DIFFERENT SPIN CASES, BECAUSE INDICES HAVE TO BE
C   FIRURED OUT FROM ITOPT2.
C
      DO 10 IRANK=1,NTOP
       CALL PCKIND(ITOPT2(IRANK),ISPIN,SYVEC1,SYVEC2,NSMSZ1,NSMSZ2,
     &             I(IRANK),J(IRANK),A(IRANK),B(IRANK))
10    CONTINUE
cjp
      if(isbwcc) then
cjp prepare spin indices according to ispin case
        if(ispin.eq.1) then
          isp1=1
          isp2=1
        else
          if(ispin.eq.2) then
            isp1=2
            isp2=2
          else
            isp1=1
            isp2=2
          endif
        endif
        write(*,202)type,SPCASE(ISPIN),iref
        else
      WRITE(*,200)TYPE,SPCASE(ISPIN)
        endif
200   FORMAT(T3,' Largest ',A1,'2 amplitudes for spin case ',A2,':')
202   format(T3,' Largest ',A1,'2 amplitudes for spin case ',A2,
     &    ' of reference no.: ',i3)
      IF(ISPIN.EQ.3)WRITE(*,299)
      WRITE(*,300)
299   FORMAT(3X,3(' ',3X,'_',3X,' ',3X,'_',13X))
300   FORMAT(3X,3('i',3X,'j',3X,'a',3X,'b',13X))
      WRITE(*,400)
400   FORMAT(77('-'))
cjp
      if(isbwcc) then
cjp ... this is not OK for uhf-cc, ignore for now
       if(bwgossip) then
      write(*,502) (irank,i(irank),j(irank),a(irank),b(irank),
     &   invhnum(i(irank),iref,isp1),
     & invhnum(j(irank),iref,isp2),invpnum(a(irank)-nocca,iref,isp1),
     & invpnum(b(irank)-nocca,iref,isp2),topt2(irank),irank=1,ntop)
        else
        write(*,500) (invhnum(i(irank),iref,isp1),invhnum(j(irank),
     &  iref,isp2),invpnum(a(irank)-nocca,iref,isp1),
     & invpnum(b(irank)-nocca,iref,isp2),topt2(irank),irank=1,ntop)
        endif
      else
      WRITE(*,500)(I(IRANK),J(IRANK),A(IRANK),B(IRANK),
     &                 TOPT2(IRANK),IRANK=1,NTOP)
      endif
      WRITE(*,400)
      WRITE(*,600)TYPE,SPCASE(ISPIN),NT2SIZ,T2NORM
      WRITE(*,400)
cpn
      if(isbwcc.and.bwgossip) then
c This will print all the amplitudes and their indexes in order
      write(*,203)type,spcase(ispin),iref
 203  format(t3,' List of all ',A1,'2 amplitudes for spin case ',A2,
     &    ' of reference no.: ',i3)
      if(ispin.eq.3)write(*,299)
      write(*,300)
      write(*,400)
      do 11 irank = 1, nt2siz
         call pckind(irank,ispin,syvec1,syvec2,nsmsz1,nsmsz2,
     &        i(1),j(1),a(1),b(1))
         itrue=invhnum(i(1),iref,isp1)
         jtrue=invhnum(j(1),iref,isp2)
cjp ... this is not OK for uhf-cc, ignore for now
         atrue=invpnum(a(1)-nocca,iref,isp1)
         btrue=invpnum(b(1)-nocca,iref,isp2)
         write(*,502)irank,i(1),j(1),a(1),b(1),itrue,jtrue,atrue,
     &         btrue,t(irank)
 11   continue
 502  format(3x,i5,6x,4i5,': true ind : ',4i5,':',f12.8)
      endif
cpn end
600   FORMAT(T3,' Norm of ',A1,'2',A2,' vector (',I9,
     &       ' symmetry allowed elements):',F14.10,'.')
500   FORMAT((2('[',I3,1X,I3,1X,I3,1X,I3,']',F8.5,1X),'[',I3,1X,
     &I3,1X,I3,1X,I3,']',F8.5))
      RETURN
      END
