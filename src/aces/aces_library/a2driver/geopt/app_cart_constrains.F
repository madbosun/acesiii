      Subroutine App_cart_constrains(HES, TMP, PMAT, GRD, GRDTMP,
     &                               CONSTRAINS, NUM_CONSTRAINS, NX)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER CONSTRAINS

      DIMENSION HES(NX, NX), TMP(NX, NX), PMAT(NX,NX),
     &          GRD(NX), GRDTMP(NX), CONSTRAINS(NUM_CONSTRAINS)
C
C Build the projected Hessian and gradients for constrained Cartesian
C optimizations. Ajith Perera, 06/2013.
C
      Call Dzero(Pmat, Nx*Nx)
      
      Do Idim = 1, Nx
            Pmat(Idim,Idim) = 1.0D0
      Enddo

      Do Idim = 1, NUM_CONSTRAINS
         Indx = Constrains(Idim) 
         Pmat(Indx, Indx) = 0.0D0
      Enddo

#ifdef _DEBUG_LVL0
      call xgemm("N", "T", Nx, Nx, Nx, 1.0D0, Pmat, Nx, Pmat, Nx,
     &            0.0d0, tmp, nx)
      write(6,"(a)") "The test for the projector: P^P"
      call output(tmp, 1, nx, 1, nx, nx, nx, 1)
      Write(6,"(a)") "The gradients, g"
      Write(6,"(5(1x,F15.7))") (Grd(I), I=1,NX)
      Write(6,"(a)") "The Hessian: H"
      call output(hes, 1, nx, 1, nx, nx, nx, 1)
#endif
C
C Project the Hessian
C
      CALL XGEMM("N", "N", NX, NX, NX, 1.0D0, PMAT, NX, HES, NX,
     &            0.0D0, TMP, NX)
      CALL XGEMM("N", "N", NX, NX, NX, 1.0D0, TMP, NX, PMAT, NX,
     &            0.0D0, HES, NX)

C Also project the gradients ie. form PG. Then copy the projected
C gradient  to the FI array so that the rest of the code can proceed.
C
      CALL XGEMM('N', 'N', NX, 1, NX, 1.0D0, PMAT, NX, GRDTMP, 
     &            NX, 0.0D0, GRD, NX)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The projected gradients: Pg"
      Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NX)
      Write(6,"(a)") "The projected Hessian: PHP"
      call output(HES, 1, nx, 1, nx,  nx, nx, 1)
#endif
      
      Return
      End
