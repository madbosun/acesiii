      SUBROUTINE ASSIGN_XYZ_LABELS (TOTNO_XYZ, XYZ_LABEL)
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Integer TOTNO_XYZ
C
C Assign a label for xyz coordiantes of each atom
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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
   
