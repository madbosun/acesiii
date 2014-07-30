      SUBROUTINE DUMPCORD(NATOMS,Q,IATNUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


      DIMENSION Q(3*NATOMS),IATNUM(NATOMS)
      ATOB = 0.529177249D0
      WRITE(6,100)
      WRITE(6,101)
      WRITE(6,100)
100   FORMAT(T2,64('-'))
101   FORMAT(T2,'Z-matrix',T13,'Atomic',T31,'C o o r d i n a t e s',
     &       /,
     &       T3,'Symbol',T13,'Number',T30,'X',T45,'Y',T60,'Z')
      IOFF=1
      DO 10 I=1,NATOMS
       LAST=linblnk(ZSYM(I))
       WRITE(6,1000)ZSYM(I)(1:LAST),IATNUM(I),(Q(J)*ATOB,
     &              J=IOFF,IOFF+2)
       IOFF=IOFF+3
10    CONTINUE
      WRITE(6,100)
1000  FORMAT(T6,A,T14,I3,T22,F14.8,T37,F14.8,T52,F14.8)
      RETURN
      END
