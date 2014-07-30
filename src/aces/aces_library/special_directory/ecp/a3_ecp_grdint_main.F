#include "flags.h"
      Subroutine A3_ECP_grdint_main(Iecp, Natoms, Namat, Iqmstr,
     &                              Jcostr, Nucstr, Nrcstr, Jstrt, 
     &                              Jrs, Expa, Nalpha, Cont, Npcoef, 
     &                              Cord, Charge, IGenby, Spherical,
     &                              Ecpgrdint, Naobasis)
C
      Implicit Double Precision(A-H, O-Z)
C
#include "ecp.par"
#include "ecp.com"
C
      Logical Spherical
      Character*2 mtype(Max_centers)
      Character*80 Cnick(Max_centers)
      Character*80 Ecpnam(Max_centers)
      Character*4 namat(max_centers)

      Dimension iqmstr(Max_centers),jcostr(Max_centers,maxjco),
     &          nucstr(Max_centers,Maxang,Maxjco),
     &          nrcstr(Max_centers,Maxang,Maxjco),
     &          Jstrt(Max_shells),Jrs(Max_shells),Expa(Nalpha),
     &          Cont(Npcoef),Cord(3,Max_centers),
     &          Charge(Max_centers), IGenby(Max_centers)

      Dimension Ecpgrdint(Natoms, 3, *)
C
CSSS      Call A3_Ecp_set4_grdint(Iecp,Nmax,namat,iqmstr,jcostr,nucstr,
CSSS     &                        nrcstr, Jstrt, Jrs, Expa, Nalpha,
CSSS     &                        Cont, Npcoef, Cord, Charge)
CSSS
CSSS      Call prepareecp(Cord, Charge, Igenby, Natoms, Ntotatoms)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The symmetry equivalent atoms"
      Write(6,"(6(1x,I3))") (IGenby(i), i=1, Ntotatoms)
      Write(6,*)
      Do j=1, Ntotatoms
      Write(6, "(3(1x,F10.6))") (Coord(i,j), i=1, 3)
      Enddo
      Write(6,*)
#endif
C
      Ntotatoms = Natoms

cVFL start 
c     Open(Unit=62, File="decpint.dat", Form="unformatted")
c     Rewind(62)

c     Write(62) Ntotatoms
cVFL end  

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The ECP deriavative integrals" 
#endif
      Call Dzero(Ecpgrdint, Natoms*3*Naobasis*(Naobasis+1)/2)

      Call A3_ECP_grdint_driver(Ecpgrdint, Ntotatoms, Natoms, 
     &                          Iloc, Isoc, IGenby, Cord, .True.,
     &                          Spherical)

cVFL start 
c     close(62)
cVFL end  

      Nint = Iloc
      If (SPherical) Nint = Isoc 
      Call write_ecpgrdint_2disk(Ntotatoms, Nint, Ecpgrdint)

      Return
      End
