      SUBROUTINE ASSIGN_XYZ_LABELS (TOTNO_XYZ, XYZ_LABEL)
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Integer TOTNO_XYZ
C
C Assign a label for xyz coordiantes of each atom
C
#include "mxatms.par"
C
      CHARACTER*4 XYZ_LABEL(TOTNO_XYZ)
      CHARACTER*3 LABEL(3*MXATMS)
C
      do i = 1, min(9,3*MXATMS)
         write(label(i),'(2i1)') 0, i
      end do
      do i = 10, min(99,MXATMS)
         write(label(i),'(i2)') i
      end do
      do i = 100, 3*MXATMS
         write(label(i),'(i3)') i
      end do
C
      DO I=1,TOTNO_XYZ
         XYZ_LABEL(I)='R'//LABEL(I)
      END DO
C            
      RETURN
      END
   
