*fordeck gauss $Revision: 1.1 $
      SubRoutine Gauss(n,lDim,A,X,C)
************************************************************************
*                                                                      *
*     purpose: Solve a set of linear equations using Gauss method      *
*                            A*X = C                                   *
*                                                                      *
*     input:                                                           *
*       A,C     : input matrices                                       *
*       n       : size of the set of equations                         *
*       lDim    : leadin dimension of A                                *
*                                                                      *
*     output:                                                          *
*       X       : solutions                                            *
*                                                                      *
*     called from: Diis, MinDns                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8 A(lDim,lDim),X(lDim),C(lDim)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*     Call qEnter('Gauss')
*
      Write(6,*) "The A matrix in Gauss", (X(I), I=1, ldim)
      Do 10 i = 1, N
         X(I) = C(I)
   10 Continue
      Do 100 i = 1, n - 1
         k = i
         Do 110 j = i + 1, n
            If (Abs(A(k,i)) .lt. Abs(A(j,i))) k = j
  110    Continue
         If (k .ne. i) Then
*           Write(*,'(A,2I3)') ' Swapping:',i,k
            Do 120 j = i, n
               Swap   = A(i,j)
               A(i,j) = A(k,j)
               A(k,j) = Swap
  120       Continue
            Swap = X(i)
            X(i) = X(k)
            X(k) = Swap
         End If
         Do 130 k = i + 1, n
            Fact = A(k,i)/A(i,i)
            Do 131 j = i + 1, n
               A(k,j) = A(k,j) - Fact*A(i,j)
  131       Continue
            X(k) = X(k) - Fact*X(i)
  130    Continue
  100 Continue
      X(n) = X(n)/A(n,n)
      Do 200 i = n - 1, 1, -1
         Do 210 k = i + 1, n
            X(i) = X(i) - A(i,k)*X(k)
  210    Continue
         X(i) = X(i)/A(i,i)
  200 Continue
*
*     Call qExit('Gauss')
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
