      Subroutine Ln_search(Energy, Cur_geo, Prv_geo_stp, Cur_grad, 
     &                     Prv_grad, Prv_geo_stpn, Nxm6, EPS)


      Implicit Double Precision (A-H, O-Z)

      Dimension Energy(2), Cur_geo(Nxm6), Prv_geo_stp(Nxm6), 
     &          Cur_grad(Nxm6), Prv_grad(Nxm6), Prv_geo_stpn(Nxm6)
 
      Data Dnull, Done, Ione, Five, Two, Three, Four, Pt1, Pt5, Pt9, 
     &            Pt15 /0.0D0, 1.0D0, 1, 5.0D0, 2.0D0, 3.0D0, 4.0D0, 
     &                 0.1D0, 0.5D0, 0.9D0, 0.15D0 /
C
      Step_Norm = Dsqrt(Ddot(Nxm6, Prv_geo_stp, 1, Prv_geo_stp, 1))

      Do Imode = 1, Nxm6
         Prv_geo_stpn(Imode) = - Prv_geo_stp(Imode)/Step_Norm
      Enddo
     
      Prv_grad_on_stp = Ddot(Nxm6, Prv_geo_stpn, 1, Prv_grad, 1)
      Cur_grad_on_stp = Ddot(Nxm6, Prv_geo_stpn, 1, Cur_grad, 1)

      Thrs = Eps*100

      If (Dabs(Prv_grad_on_stp) .lt. Thrs .or. Dabs(Cur_grad_on_stp)
     &    .lt. Thrs) Then
         Write(6, "(a,a)") "Too close to the minimum, line search",
     &          " skipped."
         Return
      Endif 
C
c The energy and gradient from the current and last point are fitted
C to a quadratic polynomial.

      Ca = Dabs(Five*Prv_grad_on_stp -Two*Cur_grad_on_stp + Three*
     &     Energy(1) - Three*Energy(2))
      Cb = Cur_grad_on_stp + Prv_grad_on_stp + Two*Energy(1) - Two*
     &     Energy(2) - Two*Ca
      Cc = Energy(2) - Energy(1) - Prv_grad_on_stp - Ca - Cb
      Cd = Prv_grad_on_stp
      Ce = Energy(1)

#ifdef _DEBUG_LVL0
      Write(6,*) "@-Ln_search, Starting parameters"
      Write(6,*) "The normalized step increment"
      Write(6,"(3F10.5)") (Prv_geo_stpn(I), I=1,Nxm6)
      Write(6,*) "Prv. and Cur. gradient nomrs"
      Write(6,"(2F10.5)") Prv_grad_on_stp, Cur_grad_on_stp
      Write(6,*) "Ca-e coefs"
      Write(6,"(5F10.5)") Ca, Cb, Cc, Cd, Ce
#endif

      If ((Prv_grad_on_stp .lt. Dnull) .and. (Cur_grad_on_stp .gt. 
     &      Dnull)) Then
c
         Crda  = Dnull
         Crdb  = DOne
         Grda  = Prv_grad_on_stp
         Grdb  = Cur_grad_on_stp
         Iycle = Ione
C
         Crdc  = Crda + (Crdb - Crda)*Max(Pt1, Min(Pt9,Dabs(Grda/
     &             Grdb)*pt5))
         Grdc  = ((four*Ca*Crdc + Three*Cb)*Crdc +Two*Cc)*Crdc + Cd
         Gcut  = 1.0D-5*Min(Dabs(Prv_grad_on_stp), Dabs(
     &               Cur_grad_on_stp))
C
#ifdef _DEBUG_LVL
      Write(6,*) "@-Ln_search, parameters at the start of search"
      Write(6,"(1x,a,4F10.5)") "Crda-b, Grda-b: ", Crda, Crdb, 
     &                          Grda, Grdb,
      Write(6,"(1x,a,3F10.5)") "Crdc, Grdc, Gcut: ", Crdc, Grdc, Gcut
#endif   

         Do While (Dabs(Grdc) .lt. Gcut) 

            Icycle = Icycle + 1
C
            If (Grdc .Gt. Dnul) Then
                Crdb = Crdc
                Grdb = Grdc
            Else
               Crda = Crdc
               Grda = Grdc
            Endif
C
            Crdc = Crda + (Crdb - Crda)*Max(Pt1, Min(Pt9,Dabs(Grda/
     &             Grdb)*pt5))
            Grdc = ((four*Ca*Crdc + Three*Cb)*Crdc +Two*Cc)*Crdc + Cd
C
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Ln_search, parameters during the ln search"
      Write(6,"(1x,a,2F10.5)") "Crdb, Grdb: ", Crdb, Grdb 
      Write(6,"(1x,a,2F10.5)") "Crda, Grda: ", Crda, Grda 
      Write(6,"(1x,a,2F10.5)") "Crdc, Grdc: ", Crdc, Grdc
#endif   
C
            If (Icycle .gt. 200) Then
               Write(6,*)
               Write(6, "(a)") "Line search failed and it is ignored!"
               Return
            Endif
C
         Enddo 
C
C If the line serch yields almost the same point as the previous point,
C do not use it.
C       
         If (Crdc .lt. pt15) Then
             Write(6, "(a)") "The line serach is skipped!"
             Return
         Else
             Do Imode = 1, Nxm6
                Cur_geo(Imode)  = Cur_geo(Imode) + 
     &                            Prv_geo_stp(Imode)*(Done - Crdc)
                Cur_grad(Imode) = Prv_grad(Imode) + (Cur_grad(imode)  
     &                            - Prv_grad(Imoded))*Crdc
             Enddo
         Endif
C 
         Energ_intpl = (((Ca*Crdc + Cb)*Crdc + Cc)*Crdc + Cd)*Crdc
     &                            + Ce

C
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Ln_search, interpolated, geo., grad. and energy"
      Write(6,*) "The geometry"
      Write(6,"(3F10.5)") (Cur_geo(i), i=1, Nxm6)
      Write(6,*) "The gradients"
      Write(6,"(3F10.5)") (Cur_grad(i), i=1, Nxm6)
      Write(6,"(a,F10.5)") "The energy", Energ_intpl
#endif
      Else
C 
         Write(6,*)
         Write(6, "(1x,a)") "Line search is skiped!"
         Write(6,*)
C
      Endif 
C
      Return
      End
 
