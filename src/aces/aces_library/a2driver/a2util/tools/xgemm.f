
cYAU - WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING
c
c    The original xgemm filtered lda, ldb, and ldc through max(1,ld?).
c This subtly changes the behavior of dgemm depending on the other arguments.
c I refused to continue that tradition. There is quite a bit of debugging
c info in this routine which should help in finding the places where the
c leading dimensions are not properly checked.

      subroutine xgemm(transa,transb,
     &                 m,n,k,alpha,a,lda,
     &                             b,ldb,
     &                       beta, c,ldc)

c ARGUMENT LIST
      character*1 transa, transb
      integer m, n, k, lda, ldb, ldc
      double precision alpha, beta
c      double precision a(lda,*), b(ldb,*), c(ldc,*)
      double precision a, b, c

c EXTERNAL FUNCTIONS

c INTERNAL VARIABLES
c Be careful here. If you build your own generic gemm, then make sure
c you use the same integer type. Otherwise, if you link against
c system-BLAS, make sure you cast the arguments correctly.




      integer   int_m, int_n, int_k, int_lda, int_ldb, int_ldc


c ----------------------------------------------------------------------

c This quick return code was taken straight from blas/dgemm. We need it
c for cases when m=lda=0 (e.g., the number of occupied orbitals is zero
c but we didn't check for it). Rather than have dgemm crash, we just
c want to return having done nothing. This also saves us from filtering
c lda, ldb, and ldc through max(1,ld?).

c BEWARE!!! xgemm will, therefore, not crash when the other arguments
c are ill-conditioned.

      if ((m.eq.0).or.
     &    (n.eq.0).or.
     &    (((alpha.eq.(0.0d0)).or.(k.eq.0)).and.(beta.eq.(1.0d0)))
     &   ) return

c ----------------------------------------------------------------------


c ----------------------------------------------------------------------

c   o recast the arguments
      int_m   = m
      int_n   = n
      int_k   = k
      int_lda = lda
      int_ldb = ldb
      int_ldc = ldc

      call dgemm(transa,transb,int_m,int_n,int_k,
     &           alpha,a,int_lda,
     &                 b,int_ldb,
     &           beta, c,int_ldc)

      return
      end

