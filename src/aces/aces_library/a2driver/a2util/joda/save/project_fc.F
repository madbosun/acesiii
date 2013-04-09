      Subroutine projec_FC(Coords, Hess, AtmMass, Grad, Work, 
     &                     Hess_project, Threshold, Nreals, 
     &                     Move_CMass, Proj_rots, Proj_grads)
C
      Implicit Double Precision (A-H, O-Z)
C
      Double Precision MO_Inertia(3,3)
      Logical Proj_rots, Proj_grads, Move_CMass
C
      Dimension Coords(3*Nreals), Grad(3*Nreals), AtmMass(Nreals), 
     &          Hess(3*Nreals, 3*Nreals), Hess_Project(3*Nreals, 
     &          3*Nreals), Work(3*Nreals, 3*Nreals), Asym_ten(3,3,3),
     &          Hess_Projects(3*Nreals, 3*Nreals), Center_Mass(3),
     &          CM(3)
 
      Data Asym_ten/ 0.0D+00,  0.0D+00,  0.0D+00,
     &               0.0D+00,  0.0D+00, -1.0D+00,
     &               0.0D+00,  1.0D+00,  0.0D+00,
     &               0.0D+00,  0.0D+00,  1.0D+00,
     &               0.0D+00,  0.0D+00,  0.0D+00,
     &              -1.0D+00,  0.0D+00,  0.0D+00,
     &               0.0D+00, -1.0D+00,  0.0D+00,
     &               1.0D+00,  0.0D+00,  0.0D+00, 
     &               0.0D+00,  0.0D+00,  0.0D+00  /
C
C Normalize the gradinet vector. 
C
      Grad_sqr = Ddot(3*Nreals, Grad, 1, Grad, 1)   
      If (Grad_sqr .GT. Threshold)  Call Normal(Grad, 3*Nreals)
C
C Move to the center of Mass coordinate system.
C
      If (Move_CMass) Then
         Call Dzero(CM, 3)
         TotMass = 0.0D0
         Do Iatoms = 1, Nreals
            Ioff = 3*(Iatoms - 1)
            TotMass = TotMass AtmMass(Iatoms)
            Do Ixyz = 1, 3
               CM(Ioff + Ixyz)   = Coords(Ioff + 1)*AtmMass(Iatoms) 
            Enddo
         Enddo
C
         Do Ixyz = 1, 3
            If (TotMass .GT. Threshold) Then
                Center_Mass(IxYz) = CM(Ixyz)/TotMass
            Else
                Write(6, "(a)") "@-Project_FC Zero total mass"
                Call Errex
            Endif
         Enddo
C         
         Do Iatoms = 1, Nreals
            Ioff  = 3*(Iatoms - 1)
            Do Ixyz = 1, 3
               Coords(Ioff + Ixyz) = Coords(Ioff + Ixyz) - 
     &                               Center_Mass(Ixyz)
            Enddo
         Enddo
      Endif
C
C compute the inverse of inertia matrix.
c 
      Call CalMOI(NReals, Coord, AtmMass, MO_Inertia)
      Call Minv(MO_Inertia, 3, 3, Work, Det, 1.0D-8, 0, 1)
C
C Mass weigh the incomming Cartesians 
C
      Ioff = 1
      TotMass = 0.0D0
      Do Iatoms = 1, Nreals
          Sqrtmass = Dsqrt(AtmMass(Iatoms))
          Ioff = Ioff + 3*(Iatom - 1)
          Do Ixyz = 1, 3
             Coords(Ioff) = Coords(Ioff)*Sqrtmass
             Ioff = Ioff + 1 
          Enddo
      Enddo    
C
C Build the Hessian projector See, Miller, Handy and Adams, JCP, 
C 72, 99, (1980).
C
      Do Iatms = 1, Natoms
         Ioff = 3*(Iatms - 1)

         Do Jatms = 1, Iatms
          Joff = 3*(Jatms - 1)
   
            Do Iz = 1, 3
               L_index = Ioff + Iz
               Jz_end = 3
               If (Iatms .EQ. Jatms)  Jz_end = Iz
C
                 Do Jz = 1, Jz_end
                    R_index = Joff + Jz
                    Sum = 0.0D0
C
                    If (Proj_rots) Then 

                        Do Ix = 1, 3
 
                           Do Iy = 1, 3
C
                              If (Asym_ten(Ix, Iy, Iz) .NE. 0.0D0) 
     &                            Then
                                   
                                   Do Jx = 1, 3

                                      Do Jy = 1, 3
                                         If (Asym_ten(Jx, Iy, Iz) .NE. 
     &                                       0.0D0) Then
                                             Sum = Sum + 
     &                                             Coords(Ioff + IY)*
     &                                             Coords(Joff + Jy)*
     &                                            MO_inertia(Ix, Jx)*
     &                                          Asym_ten(Ix, Iy, Iz)*
     &                                          Asym_ten(Jx, Jy, Jz)
                                         Endif
                                      Enddo
                                   Enddo 
                              Endif
                           Enddo
                        Enddo
                    Endif
C 
                    Hess_project(R_index, L_index) = Sum

                    If (Proj_grads) Hess_project(R_index, L_index) = 
     &                              Hess_project(R_index, L_index) + 
     &                              Grad(R_Index)*Grad(L_index)
C
                    If (Iatms .EQ. Jatms) Then
                        Hess_project(R_index, L_index) =
     &                  Hess_project(R_index, L_index) + 
     &                  Dsqrt(AtmMass(R_index)*AtmMass(R_index))/
     &                  Totmass
                    Endif
C
                 Enddo 
            Enddo 
C
         Enddo
      Enddo
C
C Build the projector (I - P)
C
      Do Jdeg = 1, 3*Nreals
          Do Ideg = 1, 3*Nreals
             If (Ideg .Eq. Jdeg)  Hess_project(Ideg, Jdeg) = 1.0D0 -
     &                            Hess_project(Ideg, Jdeg)   
             If (Ideg .GE. Jdeg .AND. Hess_project(Ideg, Jdeg) .LT. 
     &                            Threshold) Hess_project(Ideg, 
     &                                       Jdeg) = 0.0D0
              Hess_project(Jdeg, Ideg) = Hess_project(Ideg, Jdeg)
          Enddo
      Enddo
C
#ifdef _DEBUG_LVL0
      NX = 3*Nreals
      Call Eig(Hess, Work, 1, NX, 0) 
      Write(6,*)
      Write(6,"(a)") "The eigen vectors of the Hessian"
      CALL OUTPUT(Work, 1, NX, 1, Nx, Nx, Nx, 1)
      Write(6,"(a)") "The eigenvalues of the   Hessian"
      Write(6, "(4F10.5)") (Hess(I,I), I=1, NX)
#endif
C
C Project the Hessian (I - P)H(I - P)
C
      Call Xgemm("N", "N", 3*Nreals, 3*Nreals, 3*Nreals, 1.0D0, 
     &            Hess_project, 3*Nreals, Hess, 3*Nreals, 0.0D0, 
     &            Work, 3*Nreals)
C
      Call Xgemm("N", "N", 3*Nreals, 3*Nreals, 3*Nreals, 1.0D0, 
     &            Work, 3*Nreals, Hess_project, 3*Nreals, 0.0D0, 
     &            Hess, 3*Nreals)
C
#ifdef _DEBUG_LVL0
      NX = 3*Nreals
      Call Eig(Hess, Work, 1, NX, 0) 
      Write(6,*)
      Write(6,"(a)") "The eigen vectors of the projected Hessian"
      CALL OUTPUT(Work, 1, NX, 1, Nx, Nx, Nx, 1)
      Write(6,"(a)") "The eigenvalues of the  projected Hessian"
      Write(6, "(4F10.5)") (Hess(I,I), I=1, NX)
#endif
      Return
      End
     
