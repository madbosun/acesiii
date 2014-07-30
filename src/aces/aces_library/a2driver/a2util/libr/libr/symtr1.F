#include "f_types.h"
#undef SENDTO
      SUBROUTINE SYMTR1(IRREP,NUM1,NUM2,DISSIZE,A,SCR1,SCR2,ISCR)
C
C  THIS ROUTINE TRANSPOSES THE LAST TWO INDICES OF A FOUR INDEX ARRAY
C  WITHIN ONE GIVEN IRREP :
C
C   A(I,J,K,L) --> A(I,J,L,K)    WITH G(K)*G(L) = IRREP
C
C  THE TRANPOSITION IS DONE IN PLACE AND REQUIRES ONLY TWO SCRATCH
C  VECTORS OF LENGTH DISSIZE (ONE DOUBLES AND ONE INTEGERS)
C
C  INPUT :  IRREP     ...  IRREP OF G(K)*G(L)
C           NUM1      ...  POPULATION VECTOR FOR K
C           NUM2      ...  POPULATION VECTOR FOR L
C           DISSIZE   ...  DISTRIBUTION SIZE
C           A         ...  HOLDS THE MATRIX A         
C           SCR1,ISCR ...  TWO SCRATCH ARRAYS OF SIZE DISSIZE
C           SCR2      ...  (IGNORED)
C
C  OUTPUT :  A     ....  TRANSPOSED MATRIX
C
CEND
C
C  WRITTEN IN JUNE/90   JG
C  SENT THE DATA AND INDEX ARRAY TO RESORT_COLS_TO_INDEX 4/2005 ADY
C
      IMPLICIT NONE
      INTEGER IRREP,NUM1(*),NUM2(*),DISSIZE,ISCR(*)
      DOUBLE PRECISION A(DISSIZE,*),SCR1(*),SCR2(*)
#include "syminf.com"
      INTEGER iOffLK(8),iOffKL(8)
      INTEGER IRREPI,I
      INTEGER IRREPJ,J
      INTEGER IRREPK,K,NUMK
      INTEGER IRREPL,L,NUML
      INTEGER NTOTAL,IND1,IND2,IERR
      LOGICAL ICONT
      F_ADR    c_adr
      external c_adr

      if (.false.) then
         print *, '@SYMTR1: Argument analysis'
         print *, '         KL irrep ',irrep
         print *, '         K pop = ',(NUM1(K),K=1,nirrep)
         print *, '         L pop = ',(NUM2(L),L=1,nirrep)
         print *, '         lda   = ',DISSIZE
         print *, '         &A    = ',c_adr(a)
         print *, '         &SCR1 = ',c_adr(SCR1)
         print *, '         &SCR2 = ',c_adr(SCR2)
         print *, '         &ISCR = ',c_adr(ISCR)
      end if
C
      NTOTAL=0
      iOffLK(1)=0
      iOffKL(1)=0          
      DO IRREPJ=1,NIRREP-1
         IRREPI=DIRPRD(IRREPJ,IRREP)
         iOffLK(IRREPJ+1)=iOffLK(IRREPJ)+NUM2(IRREPI)*NUM1(IRREPJ)
         iOffKL(IRREPJ+1)=iOffKL(IRREPJ)+NUM1(IRREPI)*NUM2(IRREPJ)
         NTOTAL=NTOTAL+NUM1(IRREPI)*NUM2(IRREPJ)
      END DO
      IRREPJ=NIRREP
      IRREPI=DIRPRD(IRREPJ,IRREP)
      NTOTAL=NTOTAL+NUM1(IRREPI)*NUM2(IRREPJ)
      do i = 1, ntotal
         iscr(i) = i
      end do
C
C  GET FIRST THE NEW ADDRESSES AND STORE THEM IN ISCR
C
      DO 20 IRREPL=1,NIRREP
         IRREPK=DIRPRD(IRREPL,IRREP)
         NUML=NUM2(IRREPL)
         NUMK=NUM1(IRREPK)
         DO 15 L=0,NUML-1
         DO 15 K=0,NUMK-1
            IND1=1+iOffKL(IRREPL)+L*NUMK+K
            IND2=1+iOffLK(IRREPK)+K*NUML+L
            if (ind1.gt.ntotal.or.ind2.gt.ntotal) then
               print *, 'writing out of bounds: iscr(',ind1,')=',ind2
            end if
#ifdef SENDTO
            ISCR(IND1)=IND2
#else
            ISCR(IND2)=IND1
#endif
15       CONTINUE
20    CONTINUE

c      print *, 'symtr1 in:  iScr = ',(iScr(i),i=1,nTotal)
#ifdef SENDTO
      call sort_cols_to_index(dissize,ntotal,a,dissize,scr1,iscr,iErr)
#else
      call sort_cols_from_index(dissize,ntotal,a,dissize,scr1,iscr,iErr)
#endif
c      print *, 'symtr1 out: iScr = ',(iScr(i),i=1,nTotal)
      if (iErr.ne.0) then
         print *, '@SYMTR1: sort failed'
         print *, '         iErr = ',iErr
         call errex
      end if

      return
c     end subroutine symtr1
      end

