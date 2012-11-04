      Subroutine Ecp_int(Ecpint, Ecpint_4shell, La, Lb, Iloc, Jloc,
     &                   Numcoi, Numcoj)
      
      Implicit Double Precision(A-H, O-Z)
     
      Dimension Ecpint(*), Ecpint_4shell(*)

#include "ecp.par"
#include "ecp.com"

      Ind = 0
      Do Iang = 1, (La+1)*(La+2)/2
         Do Jang = 1, (Lb+1)*(Lb+2)/2
C
             Do Icon = 1, Numcoi
                Do Jcon = 1, Numcoj

                   Ind  = Ind  + 1
                   Indx = Iloc + (Icon-1)*(La+1)*(La+2)/2 + Iang
                   Jndx = Jloc + (Jcon-1)*(Lb+1)*(Lb+2)/2 + Jang

                   If (Indx .GE. Jndx) Then
                      Indecp = Indx*(Indx-1)/2 + Jndx
                      Ecpint(indecp) = Ecpint_4shell(ind)
                   Endif

                Enddo
             Enddo
C
         Enddo
      Enddo
C
      Return
      End
