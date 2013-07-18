      Subroutine Ecp_int_typ1(Iecp_cnt, Lamu, Ltot, La, Lb, Exp12, 
     &                        Xhat, Yhat, Zhat, Xp, Yp, Zp, Xc, Yc, 
     &                        Zc, Acs_xyz, Bcs_xyz, Pc2, Fact_ab, 
     &                        PC_zero, Cint, Int, Rad1_zero, Grads)

      Implicit Double Precision (A-H, O-Z)

#include "ecp.par"
#include "ecp.com"
   
      Dimension ACS_xyz(0:2*Maxang,3),BCS_xyz(0:2*Maxang,3)
      Dimension Rad(0:2*Maxang,0:2*Maxang), 
     &          Ang(0:2*maxang,0:2*Maxang),Cint(Maxmem_a3)
      
      Logical PC_zero, Rad1_zero, Grads

      Kcrl  = Kfirst(1, Iecp_cnt)
      Kcru  = Klast(1, Iecp_cnt)
C
C This is the U(lmax) term (kcrl to kcru keep track how many terms
C in the Sum Sum{j}D_{jl}r^njexp(-zeta_J)r^2) and the Clp, nlp and 
C Zlp arrays have D_{jl}, nj and Zeta_j respectively. 
C
      Call Dzero(Rad,(2*Maxang+1)**2)
C
#ifdef _DEBUG_LVL0
CSS      write(6,*) "Initialize debug values for Kcrl and Kcru"
CSS      Kcrl = 1
CSSS      Kcru = 1
      Write(6, "(a, 2(1x,i2))") "Debug values of Kcrl:Kcru : ", 
     &      Kcrl, Kcru
#endif
C
      Do Kcr = Kcrl, kcru
       
         Zeta = Zlp(Kcr)
         Np   = Nlp(Kcr)
         Dlj  = Clp(Kcr)

#ifdef _DEBUG_LVL0
      Write(6,*) "Initialize debug values for Zeta, Np and Dlj"
CSSS      Zeta = 3.0D0
CSSS      NP   = 3
CSSS      Dlj  = 2.0D0
      Write(6, "(a,(1x,i2,F10.6,1x,F10.6))") "Np,Zeta,Dlj: ",
     &      Np, Zeta, Dlj
#endif
C
         Call Ecp_rad_int_typ1(Lamu, Ltot, Np, Dlj, Exp12, Zeta, 
     &                         Xp, Yp, Zp, Xc, Yc, Zc, Pc2,
     &                         Fact_ab, Rad, Rad1_zero)

      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
#endif

      If (.Not. Rad1_zero) Call Ecp_ang_int_typ1(Xhat, Yhat, Zhat,
     &                                           La, Lb, Ltot, Lamu, 
     &                                           Acs_xyz, Bcs_xyz, 
     &                                           Rad, Ang, Cint, Int,
     &                                           Grads)

#ifdef _DEBUG_LVL0
      If (.Not. Rad1_zero) Then
      Write(6,"(a)") "The Type1 integral"
      Write(6,"(6(1x,F15.13))")(Cint(i), i=1,Int)
      Endif
#endif

      Return
      End
