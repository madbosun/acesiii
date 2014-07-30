      subroutine do_oed_to_vmol(nbasis, erd_index, scalars, ovlp_oed,
     &                          ovlp_vmol)
c-----------------------------------------------------------------------------
c   Using the shell angular momentum and basis function information, this
c   subroutine calculates an array which maps the ACES integral order to 
c   the ERD integral order.  This array may be used to re-order a block 
c   of integrals calculated by the ERD package into a format corresponding
c   to the ACES (VMOL-based) integrals.
c
c-----------------------------------------------------------------------------
      implicit none

      integer erd_index(*)
      double precision scalars(*)
      double precision ovlp_vmol(nbasis,nbasis),Ovlp_oed(nbasis,nbasis)

      integer i, j, k, ii, jj, kk  
      integer nbasis, n


      Write(6,"(6(1x,i5))") (erd_index(i), i=1, nbasis)
      Write(6,*)
      Write(6,"(6(1x,F10.5))") (scalars(i), i=1, nbasis)


      do i = 1, nbasis
         do j = 1, nbasis 

             ii = erd_index(i) 

             Ovlp_vmol(ii,j) = ovlp_oed(i,j)/scalars(ii)







      enddo
      enddo


      write(6,"(a)") "The reordered and renormalized SCF veactors"
      call output(ovlp_vmol, 1, nbasis, 1, nbasis, nbasis, nbasis,1)

      return
      end
