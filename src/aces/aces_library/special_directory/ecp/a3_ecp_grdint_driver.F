      Subroutine A3_Ecp_grdint_driver(Ecpgrdint, Ntotatoms, Natoms, 
     &                                Iloc, IGenby, Coord, Grads)

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
      Logical Grads, Zero_int

      dimension Cint(Maxmem_a3), Ecpint_4Shell(Maxints_4shell),
     &          Ecpgrd_x(Maxints_4shell), Ecpgrd_y(Maxints_4shell), 
     &          Ecpgrd_z(Maxints_4shell),
     &          IGenby(Ntotatoms), Coord(3,Ntotatoms),
     &          Ecpgrdint(Ntotatoms,3,*)
C
#ifdef _DEBUG_LVL0
CSS      Nshell = 2
      Write(6,"(a,2I3)") "The Debug value for Nshell: ", Nshell,
     &                    Ntotatoms
      do  i=1,nshell
         write(*,*) 'nshell=',i
         write(*,3913) i,katom(i),i,ktype(i)
3913     format('katom<',i2,'>  ',i2,' ktype<',i2,'>  ',i2)
         write(*,3914) i,kprim(i),i,kbfn(i)
3914     format('kprim<',i2,'>  ',i2,' kbfn<',i2,'>   ',i2)
         write(*,3915) i,nprims(i),i,ndegen(i)
3915     format('nprims<',i2,'> ',i2,' ndegen<',i2,'> ',i2)
         write(*,3916) i,kmini(i),i,kmaxi(i)
3916     format('kmini<',i2,'>  ',i2,' kmaxi<',i2,'>  ',i2)
         write(*,*)
      enddo
#endif
C
      Iloc   = 0
      Do Icent = 1, Ntotatoms 
      do ishell = 1, Nshell
CSSS      do ishell = 1, 1
C
        iatom  = katom(ishell)
        If (IGenby(Icent) .EQ. Iatom) Then

        xa = Coord(1,Icent)
        ya = Coord(2,Icent)
        za = Coord(3,Icent)
        la = ktype(ishell) - 1

        kprimi=kprim(ishell)-1
        idegen=ndegen(ishell)
        idegen_grd=Ideg_grd(la)
        numcoi=numcon(ishell)

C begin loop jshell 
 
        Jloc = 0
        Do Jcent = 1, Ntotatoms
        do jshell=1, Nshell
CSSS        do jshell=1, 1

C
C get symmetry infor regarding shell j
C
          jatom=katom(jshell)
          lb = ktype(jshell) - 1
          Jdegen=ndegen(jshell)
C
          If (IGenby(Jcent) .EQ. Jatom) Then
 
          xb=Coord(1,Jcent)
          yb=Coord(2,JCent)
          zb=Coord(3,Jcent)

          kprimj=kprim(jshell)-1
          numcoj=numcon(jshell)

C Get the density for this shell pair.
     
CSSS          Call Dzero(Dens_4shell, Maxints_4shell)
CSSS          Call Get_shell_den(Dens_fao, Dens_4shell, Naobfns, La, Lb,
CSSS     &                       Iloc, Jloc, Numcoi, Numcoj)

          Do iprim = 1, nprims(ishell)
CSSS          Do iprim =1, 1

             exp1 = expnt(kprim(ishell)+iprim-1)
             Indnpri = Kprimi + Iprim
                     
            Do jprim = 1, nprims(jshell)
CSSS            Do jprim = 1, 1 

CSS                Call Dzero(Ecpint_4shell, Maxints_4shell)
CSS                Call Dzero(Cint, Maxmem_a3)

                exp2 = expnt(kprim(jshell)+jprim-1)
                Indnprj = Kprimj + Jprim
C
#ifdef _DEBUG_LVL0
      Write(6, "(a, 4(1x,I2))") "The Ishell, Jshell, La, Lb: ", 
     &                           Ishell, Jshell, La, Lb
      Write(6, "(a, 4(1x,I2))") "The prim. pair; Iprim, Jprim",
     &                           Iprim, Jprim
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
                Do Icnt = 1,  Ntotatoms
                   Isym_unq_cnt = IGenby(Icnt)
                   If (Ipseux(Isym_unq_cnt) .NE. 0) Then
                      
                      If (.NOT. (Grads .AND. Icent .EQ. Icnt)) Then
                         Iecp_cnt = Isym_unq_cnt
                         Call Dzero(Ecpint_4shell, Maxints_4shell)
                         Call Dzero(Cint, Maxmem_a3)
 
                         Call Ecp_grdint_4prim(Xa, Ya, Za, Xb, Yb, 
     &                                         Zb, La, Lb, Coord, Exp1,
     &                                         Exp2, Natoms, Ntotatoms, 
     &                                         Cint, Int, Zero_int, 
     &                                         Icnt, Iecp_cnt, Grads)

                         If (.NOT. Zero_int)
     &                   Call ecp_int_4shell(Cint, Ecpint_4shell, 
     &                                       La, Lb, Numcoi, Numcoj, 
     &                                       Indnpri, Indnprj, Iprim,
     &                                       Jprim, Jnt, Grads)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,6(1x,I2))")
     &"The Contracted integral for shell pair", Ishell, JShell, Iprim,
     & Jprim, Idegen_grd, Jdegen
       Do I = 1, Idegen_grd
          Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,Jdegen)
      Enddo
#endif

                         Call A3_Ecpgrd_int(Ecpint_4shell, La, Lb, 
     &                                      Idegen_grd, Idegen,
     &                                      Jdegen, Numcoi, Numcoj,
     &                                      Iprim, Jprim, Exp1, Exp2,
     &                                      Ecpgrd_x, Ecpgrd_y, 
     &                                      Ecpgrd_z, Ntotatoms, 
     &                                      Icent, ICnt)
                        
                         Call Ecp_grdint(Ecpgrdint, Ecpgrd_X, Ecpgrd_Y,
     &                                   Ecpgrd_Z, La, Lb, Iloc, Jloc, 
     &                                   Numcoi, Numcoj, Icent, Icnt,
     &                                   Ntotatoms)

                      Endif
                   Endif 
C
                Enddo
C
             Enddo
          Enddo        

          Jloc = Jloc + (Lb+1)*(Lb+2)*Numcoj/2

        Endif 
        Enddo
        Enddo
        Iloc = Iloc + (La+1)*(La+2)*Numcoi/2

      Endif 
      Enddo 
      Enddo

      return
      end
