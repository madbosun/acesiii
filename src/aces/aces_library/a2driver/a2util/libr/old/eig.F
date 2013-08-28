      SUBROUTINE EIG(A,B,L,N,N1)
C
C Shell to drive Householder algorithm.
C A - MATRIX TO BE DIAGONALIZED (EIGENVALUES IN DIAGS AFTERWARDS)
C B - EIGENVECTORS RETURNED IN COLUMNS
C L - not used
C N - SIZE OF MATRIX
C N1 - EIGENVECTORS AND EIGENVALUES ARE REORDERED, with eigenvalues:
c      1     -- unordered (currently ascending because eispack does it)
c      0     -- ascending
c      other -- descending
C 
CEND
c 91/10/1 DCC use eispack routines instead of numerical recipes
C
      integer l, n, n1
      double precision A(N,N),B(N,N)
c
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
c
      integer i, j, mu, ierr
      double precision w1, w2
C
      IF(N.EQ.1)THEN
         B(1,1)=1.0
         RETURN
      ENDIF
c
      CALL DCOPY(N*N,A,1,B,1)
      CALL ZERO(A,N*N)
      CALL TRED2(N,N,B,A(1,1),A(1,2), b)
      CALL TQL2(N,N,A(1,1),A(1,2),B, ierr )
      if ( ierr .ne. 0 ) then
         write (6,*) 'eig:  eigenvalue not found -- ', ierr
         call errex
      endif
c
      DO 103 I=1,N
         A(I,I)=A(I,1)
 103  CONTINUE
      DO 301 I=1,2
      DO 301 J=I+1,N
         A(J,I)=0.D0
 301  CONTINUE
      a(1,2) = 0.0d0
c
c make first significant element of each eigenvector positive
c
      tol=1.d-5
      do 31 i=1,n
       jfirst=0
       do 32 j=1,n
        if(abs(b(j,i)).gt.tol.and.jfirst.eq.0)jfirst=j 
32     continue
       if(b(jfirst,i).lt.0.d0)call vminus(b(1,i),n)
31    continue
 
c
      IF(N1.EQ.1 .or. n1.eq.0 ) RETURN
c
      DO 20 I=1, n-1
      DO 20 J=I+1,N
         IF ( A(I,I) .lt. A(J,J) ) then
            W1=A(I,I)
            A(I,I)=A(J,J)
            A(J,J)=W1
            DO 10 MU=1,N
               W2=B(MU,I)
               B(MU,I)=B(MU,J)
               B(MU,J)=W2
 10         continue
         endif
 20   CONTINUE
      RETURN
      END
