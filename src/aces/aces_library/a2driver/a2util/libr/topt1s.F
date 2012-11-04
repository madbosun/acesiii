      SUBROUTINE TOPT1S(ISPIN,ILIST,NTOP,NT1SIZ,NOCC,
     &                  NVRT,NSMSZ,T,TOPT1,
     &                  SYVEC,ITOPT1,I,A,TYPE)
C
C THIS ROUTINE PICKS UP A T1 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT1 AND THEIR ASSOCIATED INDICES
C  IN THE I AND A VECTORS.
C
C INPUT:
C       ISPIN - THE LIST SUBTYPE
C       ILIST - THE LIST NUMBER
C       NTOP  - THE NUMBER OF LARGEST AMPLITUDES (BY ABSOLUTE VALUE)
C               WHICH ARE TO BE FOUND
C       NT1SIZ- THE TOTAL SIZE OF THE T1 VECTOR
C       NOCC  - THE NUMBER OF OCCUPIED ORBITALS OF THIS SPIN TYPE
C       NVRT  - THE NUMBER OF VIRTUAL ORBITALS OF THIS SPIN TYPE
C       NSMSZ - THE TOTAL SIZE OF THE SYMMETRY VECTOR (NOCC*NVRT)
C       SYVEC - THE SYMMETRY VECTOR FOR THE A,I OR a,i LIST
C       TYPE  - ???????????
C
C OUTPUT:
C       TOPT1 - THE NTOP LARGEST ELEMENTS IN THE T1 VECTOR
C               SORTED BY ABSOLUTE VALUE
C       I     - THE I INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       A     - THE A INDICES CORRESPONDING TO THE VALUES IN TOPT1
C
C SCRATCH:
C
C       T     - USED TO HOLD THE SYMMETRY-PACKED T VECTOR
C       ITOPT1- USED TO HOLD THE OFFSETS CORRESPONDING TO THE
C               ELEMENTS IN TOPT1
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A(NTOP),SYVEC(NSMSZ)
      CHARACTER*2 SPCASE(3)
      CHARACTER*1 TYPE
      DIMENSION T(NT1SIZ),TOPT1(NTOP),ITOPT1(NTOP),I(NTOP)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
cjp
#include "bwcc.com"
      DATA SPCASE /'AA','BB','AB'/
C
C THIS ROUTINE PICKS UP A T1 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT1 AND THEIR ASSOCIATED INDICES
C  IN THE I AND A VECTORS.
C
C
C STATEMENT FUNCTIONS
C
      IGETI(IX,NUMA)=1+(IX-1)/NUMA
      IGETA(IX,NUMA)=IX-(IGETI(IX,NUMA)-1)*NUMA
      NNM1O2(IX)=(IX*(IX-1))/2
C
      IF(NT1SIZ.EQ.0)RETURN
C
      CALL GETLST(T,1,1,1,ISPIN,ILIST)
C
C COMPUTE NORM OF T1 VECTOR
C
      X=SNRM2(NT1SIZ,T,1)
      T1NORM=X
C
C PUT NTOP LARGEST ELEMENTS INTO TOPT1 AND THEIR ASSOCIATED INDICES INTO
C   ITOPT1.
C
      CALL SCANVC(T,TOPT1,ITOPT1,NTOP,NT1SIZ)
      DO 10 IRANK=1,NTOP
       IUNPK=SYVEC(ITOPT1(IRANK))
       A(IRANK)=IGETA(IUNPK,NVRT)+NOCCO(ISPIN)
       I(IRANK)=IGETI(IUNPK,NVRT)
10    CONTINUE
cjp
      if(isbwcc) then
        write(*,201)type,SPCASE(ISPIN),iref
        else
      WRITE(*,200)TYPE,SPCASE(ISPIN)
        endif
200   FORMAT(T3,' Largest ',A1,'1 amplitudes for spin case ',A2,':')
201   format(T3,' Largest ',A1,'1 amplitudes for spin case ',A2,
     +   ' of reference no.: ',i3)
      IF(ISPIN.EQ.3)WRITE(*,299)
      WRITE(*,300)
299   FORMAT(3X,3(' ',3X,'_',3X,' ',3X,'_',13X))
300   FORMAT(3X,3('i',3X,' ',3X,'a',3X,' ',13X))
      WRITE(*,400)
400   FORMAT(77('-'))
      if(isbwcc) then
      if(bwgossip) then
        write(*,501) (irank,i(irank),a(irank),
     &         invhnum(i(irank),iref,ispin),
     &        invpnum(a(irank)-nocco(ispin),iref,ispin),
     &         topt1(irank),irank=1,ntop)
        else
         write(*,500)(invhnum(i(irank),iref,ispin),
     &      invpnum(a(irank)-nocco(ispin),iref,ispin),
     &              topt1(irank),irank=1,ntop)
        endif
      else
      WRITE(*,500)(I(IRANK),A(IRANK),
     &                 TOPT1(IRANK),IRANK=1,NTOP)
      endif
      WRITE(*,400)
      WRITE(*,600)TYPE,SPCASE(ISPIN),NT1SIZ,T1NORM
      WRITE(*,400)
600   FORMAT(T3,' Norm of ',A1,'1',A2,' vector (',I9,
     &        ' symmetry allowed elements):',F14.10,'.')
500   FORMAT((2('[',I3,1X,3X,1X,I3,1X,3X,']',F8.5,1X),'[',I3,1X,
     &3X,1X,I3,1X,3X,']',F8.5))
cjp
501    format(i4,':',i4,1x,i4,':true indices:',i4,1x,i4,' : ',f12.8)
      RETURN
      END
