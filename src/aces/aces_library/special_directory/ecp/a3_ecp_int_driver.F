      subroutine A3_Ecp_int_driver(Ecpint, ECpint_C, Natoms, Ntotatoms,
     &                             IGenby, Coord, Grads, Spherical)
c----------------------------------------------------------------------
c     Ecp_int_driver drives the calculation of 1e-matrices in the basis of
c     cartesian atomic orbitals using subroutine ecp_int.
c----------------------------------------------------------------------
c     tol    = integral cutoff (ints to be neglected, if value expected
c              to be < exp(-tol)(i)
c     Cint   = 1e-integral array ('batch') of dimension nfij
c       
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C-----------------------------------------------------------------------
#include "ecp.par"
#include "ecp.com"

C Ecpints need to be (2*Maxang+1)**2 times the square of the 
C maximum number of contracted functions per shell (currently
C set at 20 for ECP integrals). The current Max_cbf is 1000. 
C
      Logical Grads, Zero_int, Spherical

      dimension Cint(Maxmem_a3), Ecpint_4Shell(Maxints_4shell),
     &          IGenby(Ntotatoms),
     &          Coord(3,Ntotatoms),Ecpint(Max_cbf*(Max_cbf+1)/2)

c-----------------------------------------------------------------------
#ifdef _DEBUG_LVL0
CSS      Nshell = 2
      Write(6,"(a,I3)") "The Debug value for Nshell: ", Nshell
CSS      nprims(2) = 1
#endif

      Iloc   = 0
      Isoc   = 0
      Do Icent = 1, Ntotatoms
      do ishell = 1, Nshell
CSSS      do ishell = 1, 1
        iatom  = katom(ishell)

        If (IGenby(Icent) .EQ. Iatom) Then

        xa = Coord(1,icent)
        ya = Coord(2,icent)
        za = Coord(3,icent)
        la = ktype(ishell) - 1

        kprimi=kprim(ishell)-1
        idegen=ndegen(ishell)
        numcoi=numcon(ishell)

C begin loop jshell 

        Jloc = 0
        Jsoc = 0
        Do Jcent = 1, Icent
        do jshell=1, Nshell
CSSS        do jshell=1, 1

          Jdegen=ndegen(Jshell)
          jatom=katom(jshell)
          lb = ktype(jshell) - 1
C
#ifdef _DEBUG_LVL0       
      Write(6,"(a)") "Icent, Jcent, Jatom, Ishell, Jshell: " 
      write(6,"(6(1x,i4))") Icent, Jcent, Jatom, Ishell, Jshell, Jloc
      Write(6,*)
#endif
          If (IGenby(Jcent) .EQ. Jatom) Then
          If (.NOT. (ICent .EQ. Jcent .AND. JShell .GT. IShell)) Then
C
          xb=Coord(1,jcent)
          yb=Coord(2,jcent)
          zb=Coord(3,jcent)

          kprimj=kprim(jshell)-1
          numcoj=numcon(jshell)

#ifdef _DEBUG_LVLM
      Nprims(1) = 1
      Write(6,"(a,I3)") "The Debug value for Nprims: ", Nprims(1)
#endif
          Call Dzero(Ecpint_4shell, Maxints_4shell)
C
          Do iprim = 1, nprims(ishell)
CSSS                  Do iprim = 3, 3

              exp1 = expnt(kprim(ishell)+iprim-1)
              Indnpri = Kprimi + Iprim
                     
              Do jprim = 1, nprims(jshell)
CSSS                     Do jprim = 3, 3
                 exp2 = expnt(kprim(jshell)+jprim-1)
                 Indnprj = Kprimj + Jprim
                 Call Dzero(Cint, Maxmem_a3)
C
#ifdef _DEBUG_LVL0
      Write(6, "(a, 4(1x,I2))") "The Ishell, Jshell, La, Lb: ", 
     &                           Ishell, Jshell, La, Lb
      Write(6, "(a, 4(1x,I2))") "The prim. pair; Iprim, Jprim",
     &                           Iprim, Jprim
#endif
#ifdef _DEBUG_LVL0
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of A : ", Xa, ya, Za
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of B : ", Xb, yb, Zb
      Write(6, "(a, 2(1x,F10.6))") "The Primtive pair: ", Exp1, Exp2
      Write(6,*)
      Write(6,*) "Initilize debug data"
#endif 
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Initilize debug xyz of A and B"
      xa = 1.0
      ya = 1.0
      za = 1.0
      xb = -1.0
      yb = -1.0
      zb = -1.0
      exp1 = 1.5
      exp2 = 1.5 
      La = 2
      Lb = 2
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of A : ", Xa, ya, Za
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of B : ", Xb, yb, Zb
      Write(6, "(a, 2(1x,F10.6))") "The Primtive pair: ", Exp1, Exp2
      Write(6,*)
#endif
                 Call Ecp_int_4prim(Xa, Ya, Za, Xb, Yb, Zb, 
     &                              La, Lb, Coord, Exp1, Exp2,
     &                              Natoms, Ntotatoms, IGenby,
     &                              Cint, Int, Zero_int, 
     &                              Grads)

                 If (.NOT. Zero_int)
     &              Call ecp_int_4shell(Cint, Ecpint_4shell, 
     &                                  La, Lb, Numcoi, Numcoj, 
     &                                  Indnpri, Indnprj, Iprim,
     &                                  Jprim, Jnt)
             Enddo
          Enddo        
C 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,4(1x,I2))")
     &"The Contracted integral for shell pair", Ishell, JShell, Idegen,
     & Jdegen

      If (Ishell .NE. Jshell) Then
          Do I = 1, Idegen
             Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,Jdegen)
          Enddo
      Else
          Do I = 1, Idegen
             Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,I)
          Enddo
      Endif
#endif  
                  Call Ecp_int(Ecpint, Ecpint_4shell, 
     &                         La, Lb, Iloc, Jloc, Isoc, Jsoc, 
     &                         Numcoi, Numcoj, Ishell, Jshell,
     &                         Spherical)
C
          Jloc = Jloc + (Lb+1)*(Lb+2)*Numcoj/2
          Jsoc = Jsoc + (2*Lb+1)*NUmcoj
C
        Endif
        Endif
C
        Enddo
        Enddo

        Iloc = Iloc + (La+1)*(La+2)*Numcoi/2
        Isoc = Isoc + (2*La+1)*Numcoi
C
      Endif 
      Enddo 
      Enddo
C     
      Ndim = Iloc
      If (Spherical) Ndim = Isoc 
#ifdef _DEBUG_LVL0
      w=0.0d0
      Write(6,"(a)") "@-Ecpint_Driver, The ECP integrals"
      Write(6,"(a,1x,i3)") "The dimension of ints: ", Ndim
      Write(6, "(6(1x,F10.7))") (Ecpint(I), I=1, Ndim*(Ndim+1)/2)
      do i =1, Ndim*(Ndim+1)/2
      W = w+ecpint(i)*ecpint(i)
      enddo
      Write(6, "(a,(1x,F20.13))") "The ecpint. check sum = ", w
#endif
C
      Open(unit=61, file='ecpint.dat', form='unformatted')
      Rewind 61
      write(61) Ndim
      Write(61) (Ecpint(i), i = 1, Ndim*(Ndim+1)/2)
      Close (61)
      Call Putrec(20, "JOBARC", "ECP1INTS", Ndim*(Ndim+1)/2, Ecpint)
C
      return
      end
