
cYAU - WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING
c
c    For some reason, which I am not privy to, the original xgeev was not bound
c the same way dgeev is.

      subroutine mn_geev(ndim,  nlow,  amat,
     &                   evalr, evali, evec,
     &                   scr,   nscr,  ierr)
      implicit none

      double precision amat(*), evalr(*), evali(*), evec(*), scr(*)
      integer ndim, nlow, ierr, nscr

      call xgeev('N', 'V', ndim, amat, nlow,
     &           evalr, evali, evec, ndim,
     &           evec, ndim, scr, nscr, ierr)

      return
      end

