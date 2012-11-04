      Subroutine Ecp_grdint(Ecpgrdint, Ecpgrd_X, Ecpgrd_Y, Ecpgrd_Z, 
     &                      La, Lb, Iloc, Jloc, Numcoi, Numcoj, 
     &                      Iatom, Jatom, Natoms)
      
      Implicit Double Precision(A-H, O-Z)
    
      Dimension Ecpgrdint(Natoms,3,*), Ecpgrd_X(*), Ecpgrd_Y(*),
     &          Ecpgrd_Z(*)

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

                   If (Indx .GE. Jndx) Indecp = Indx*(Indx-1)/2 + Jndx
                   If (Jndx .GT. Indx) Indecp = Jndx*(Jndx-1)/2 + Indx
#ifdef _DEBUGLVL0
      Write(6,"(2(1x,I3),3(1x,F10.5))") indecp, ind,
     &                               ecpgrd_x(ind), ecpgrd_y(ind),
     &                               ecpgrd_z(ind)
#endif

                   Ecpgrdint(Iatom,1,Indecp) = Ecpgrd_X(ind) +
     &                                         Ecpgrdint(Iatom,1,Indecp)
                   Ecpgrdint(Iatom,2,Indecp) = Ecpgrd_Y(ind) +
     &                                         Ecpgrdint(Iatom,2,Indecp)
                   Ecpgrdint(Iatom,3,Indecp) = Ecpgrd_Z(ind) +
     &                                         Ecpgrdint(Iatom,3,Indecp)

                   Ecpgrdint(Jatom,1,Indecp) = -Ecpgrd_X(ind) +
     &                                         Ecpgrdint(Jatom,1,Indecp)
                   Ecpgrdint(Jatom,2,Indecp) = -Ecpgrd_Y(ind) + 
     &                                         Ecpgrdint(Jatom,2,Indecp)
                   Ecpgrdint(Jatom,3,Indecp) = -Ecpgrd_Z(ind) +
     &                                         Ecpgrdint(Jatom,3,Indecp)

                Enddo
             Enddo
C
         Enddo
      Enddo
C
      Return
      End
