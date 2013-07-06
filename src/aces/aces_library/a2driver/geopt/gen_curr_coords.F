      Subroutine Gen_curr_coords(React, Prodt, Currn, Ndim, XYZ)

      Implicit Double Precision (A-H,O-Z)
 
      Logical XYZ

#include "mxatms.par"
#include "machsp.com"
#include "coord.com"

      Dimension React(Ndim), Prodt(Ndim), Currn(Ndim)

      Call DGetrec(20, "JOBARC", "RXSTRUCT", Ndim, React)
      Call DGetrec(20, "JOBARC", "PRSTRUCT", Ndim, Prodt)

CSSS      Call XDcopy(Ndim, React, 1, Currn, 1)

      If (XYZ) Then
      
         Do Iatm = 1, Ndim/3
            Do jatm = 1, Iatm
               Distij = Dist(React(3*(Iatm-1)+1), React(3*(Jatm-1)+1))
               Distkl = Dist(Prodt(3*(Iatm-1)+1), Prodt(3*(Jatm-1)+1))
               
               Distmn = Distij + (Distkl - Distij)/2.0D0
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Distij, Distkl, Distmn"
      Write(6,"(3(1x,F12.6))") Distij, Distkl, Distmn*0.529177249D0
#endif
            Enddo
         Enddo
C
         Call DPutrec(1,'JOBARC','CUSTRUCT',Ndim,Currn)
      
      Else
  
         Do idim = 1, Ndim
 
CSSS            Currn(Idim) = React(Idim) + (Prodt(Idim) - 
CSSS     &                                   React(Idim))/2.0D0
         Enddo
 
CSSS         Call XDcopy(Ndim, Currn, 1, R, 1)
         Call DPutrec(20, "JOBARC", "RXSTRUCT", Ndim, Currn)

      Endif
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The current structue for LST optimizations"
      Write(6,"(3(1x,F12.6))") (currn(i),i=1, Ndim)
#endif 
      Return
      End
