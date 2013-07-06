      Subroutine Read_ric_constrains(Constrains, 
     &                               Num_bond_Constrains, 
     &                               Num_angle_constrains,
     &                               Num_dihs_angle_constrains,
     &                               Not_found)
      Implicit double precision (a-h,o-z)
      Logical ZMAT_PRESENT
      Character*80 wrk
      Integer Constrains 
      Logical Not_found
C
#include "mxatms.par"

      Dimension Constrains(4, Maxredunco)
  
      Not_found = .False.

      INQUIRE(FILE='ZMAT',EXIST=ZMAT_PRESENT)
C
      IF (ZMAT_PRESENT)THEN
         OPEN(UNIT=4,FILE='ZMAT',FORM='FORMATTED',STATUS='OLD')
         rewind(4)
      else
         write(6,'(T3,a)')'@Read_constrains, ZMAT file not present.'
         call errex
      ENDIF
c
300   READ(4,'(A)', END=900) WRK
      IF (WRK(1:11) .NE.'*CONSTRAINS') goto 300

      READ(4,'(A)', END=900) WRK
      IF (WRK(1:11) .EQ.'*END') GOTO 999
      BACKSPACE(4)

      READ(4,*,END=900) Num_bond_Constrains
      Do Iread = 1, Num_bond_constrains
          READ(4,*,END=900) Constrains(1,Iread), Constrains(2,Iread)
      Enddo 
C
      Ibegin = Num_bond_constrains + 1
   
      READ(4,'(A)', END=900) WRK
      IF (WRK(1:11) .EQ.'*END') THEN
         Num_angle_constrains      = 0 
         Num_dihs_angle_constrains = 0
         GOTO 99
      ENDIF
      BACKSPACE(4)

      READ(4,*,END=900) Num_angle_constrains
      Do Iread = Ibegin, Ibegin+Num_angle_constrains-1
          READ(4,*,END=900) Constrains(1,Iread), Constrains(2,Iread),
     &                      Constrains(3,Iread)
      Enddo 

      Ibegin = Num_bond_constrains + Num_angle_constrains + 1

      READ(4,'(A)', END=900) WRK
      IF (WRK(1:11) .EQ.'*END') THEN
         Num_dihs_angle_constrains = 0
         GOTO 99
      ENDIF
      BACKSPACE(4)

      READ(4,*,END=900) Num_dihs_angle_constrains 
      Do Iread = Ibegin, Ibegin+Num_dihs_angle_constrains-1
          READ(4,*,END=900) Constrains(1,Iread), Constrains(2,Iread),
     &                      Constrains(3,Iread), Constrains(4,Iread)
      Enddo 

      GO TO 99
C
900   WRITE(6,901)
901   FORMAT(3x,'@Read_constrains, *CONSTRAINS namelist not found or',
     &       ' incomplete.')
      GO TO 999
C
 999  CONTINUE
      CLOSE(UNIT=4,STATUS='KEEP')
      Write(6,*)
      Write(6,"(3x,2a)") "Constrained redundent internal ",
     &                   "optimization is requested but no constraints"
      Write(6,"(3x,a)")  "are specified in the inputfile (ZMAT)."
      Write(6,"(3x,2a)") "Calculation proceeds as a unconstrained ",
     &                   "optimization."
C 
      Not_found = .True.
      Num_bond_constrains       = 0
      Num_angle_constrains      = 0
      Num_dihs_angle_constrains = 0
C
  99  CONTINUE
      CLOSE(UNIT=4,STATUS='KEEP')

      Return
      End
