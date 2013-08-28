      SUBROUTINE EIG2(A,B,L,N,N1)
C
C Shell to drive Householder algorithm.  Slightly different from 
C EIG routine in that less memory is used.
C
C B - MATRIX TO BE DIAGONALIZED (EIGENVECTORS RETURNED)
C A - N by 2 SCRATCH VECTOR.  A(I,1) CONTAINS Ith EIGENVALUE ON RETURN
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
      double precision A(N,2),B(N,N)
c
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
c
      integer i, j, mu, ierr
      double precision w1, w2
C
      IF(N.EQ.1)THEN
         A(1,1)=B(1,1)
         B(1,1)=1.0
         RETURN
      ENDIF
c
      CALL TRED2(N,N,B,A(1,1),A(1,2), b)
      CALL TQL2(N,N,A(1,1),A(1,2),B, ierr )
      if ( ierr .ne. 0 ) then
         write (6,*) 'eig:  eigenvalue not found -- ', ierr
         call errex
      endif
c
      IF(N1.EQ.1 .or. n1.eq.0 ) RETURN
c
      DO 20 I=1, n-1
      DO 20 J=I+1,N
         IF ( A(I,1) .lt. A(J,1) ) then
            W1=A(I,1)
            A(I,1)=A(J,1)
            A(J,1)=W1
            DO 10 MU=1,N
               W2=B(MU,I)
               B(MU,I)=B(MU,J)
               B(MU,J)=W2
 10         continue
         endif
 20   CONTINUE
c
      RETURN
      END
