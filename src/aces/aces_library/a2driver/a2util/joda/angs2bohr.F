      Subroutine Angs2bohr(R, NX)
C 
      Implicit Double Precision (A-H, O-Z)
      Dimension R(NX)
C
      ATOR = DACOS(-1.0D0)/180.0D0
      ATOB = 0.529177249D0
c
      Do IDeg = 8, NX-1, 3
         If (Ideg .NE. 8) R(Ideg + 1) = R(Ideg + 1)*ATOR
                          R(ideg) = R(Ideg)*ATOR
      Enddo

      Do IDeg = 4, NX-2, 3
         R(Ideg) =  R(Ideg)/ATOB
      Enddo

      Return
      End
