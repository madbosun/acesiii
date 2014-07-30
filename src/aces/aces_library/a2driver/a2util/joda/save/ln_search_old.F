
      Subroutine Ln_Search(Energy, q, g, dq, nInter, nIter, dqHdq,
     &                     Thr)
      Implicit Double Precision (a-h,o-z)

      DIMENSION Energy(2), A(0:5), Projg(2), q(nInter, 2),
     &          g(nInter, 2), dq(nInter, 2), B(4,0:3), FVal(4)
      Logical RC
      DATA Zero /0.0D0/, One /1.0D0/, Two /2.0D0/, Three /3.0D0/
     &     Four /4.0D0/, Half /0.5D0/
*
      A(0)=Zero
*
      iNew=nIter-1
      iOld=nIter
*
*     1) Quartic search, -0.5 <= Xmin <= 0.5 
*     2) Cubic search, -1.5 <= Xmin <= 1.5
*     3) E(iNew) is the best, no action
*     4) E(iNew) is not the best,
*        XMin=0.0 between iNew and the Best point
*
*---- Quartic line search
*
c      Write(6,*) "The distance array"
c      call output(q, 1, ninter, 1, 2, ninter, 2, 1)
c      Write(6,*) "The gradient array"
c      call output(g, 1, ninter, 1, 2, ninter, 2, 1)
c      Write(6,*) "the energies:", energy(1), energy(2)
      call dscal(nInter*2, -One, g, 1)
      Projg(1)=-DDot(nInter,q(1,iNew),1,g(1,iOld),1)
     &         +DDot(nInter,q(1,iOld),1,g(1,iOld),1)
      Projg(2)=-DDot(nInter,q(1,iNew),1,g(1,iNew),1)
     &         +DDot(nInter,q(1,iOld),1,g(1,iNew),1)
      A(3)=-(Two*(Energy(iNew)-Energy(iOld))
     &    +(Projg(2)+Projg(1)))
      A(1)= (Energy(iNew)-Energy(iOld))
     &    - (One/Four) * A(3)
C       Edif = -(Energy(iNew)-Energy(iOld))
C       Gmid = (Projg(1) + Projg(2))/Two
C       A(4) = -(Three*Edif + gmid)/Two
C       A(2) = Two*(Edif + Gmid)
*     Write (*,*) ' Projg(1),Projg(2),A(3),A(1)=',
*    &              Projg(1),Projg(2),A(3),A(1)
*
*---- Test the gradient difference
*
      Test1=Projg(2)-Projg(1)
*     Write (*,*) ' Test1=',Test1
      If (Test1.lt.Thr) Then
         Write (*,*) '-- Line search failed, negative curvature'
         Return
      End If
*
*---- Now check if data is consistent with the condition that
*     the second derivative is zero at the minimum.
*
      Test2=Test1**2-Three*A(3)**2
C      Test2 = Test1**2 - Three*A(2)**2 
*     Write (*,*) ' Test2=',Test2
      If (Test2.lt.Zero) Then
         Write (*,*) 
     &       '-- Quartic line search failed, nonzero 2nd derivative'
         Go To 100
      End If
*
      A(2)=(Sqrt(Test2)+Test1)/Four
C      A(3) = (DSqrt(Test2) + Test1)/Four
      If (A(2).lt.Thr) Then
         Write (*,*) 
     &       '-- Quartic line search failed, A(3) too small,',
     &       'using cubic minimum'
         Go To 100
      End If
*
      A(4)=Half*(Test1-Test2)
      nOrder=4
      XStart=Zero
      XLow=-Half
      XHi = Half
      Go To 500
*
*---- Cubic line search
*
*     Range: -0.5 - 0.5
*
 100  Continue
      Write (*,*) '-- Cubic line search'
      FVal(1)=Energy(iOld)
      FVal(2)=Energy(iNew)
      FVal(3)=Projg(1)    
      FVal(4)=Projg(2)    
*     B=Zero
      x0=One
      x1=Half
      x2=x1*Half
      x3=x2*Half
      x0g=Zero
      x1g=One
      x2g=Two*x1
      x3g=Three*x2
      Fact=One 
      ii=1
      jj=3
      Do i = 1, 2
         Fact=-One*Fact
         B(ii,0)=x0
         B(ii,1)=Fact*x1
         B(ii,2)=x2
         B(ii,3)=Fact*x3
         B(jj,0)=x0g
         B(jj,1)=x1g
         B(jj,2)=Fact*x2g
         B(jj,3)=x3g
         ii = ii + 1
         jj = jj + 1
      End Do
      nOrder=3
c      Write(6,*) "Before The call to Gauss", (FVAL(I), I=1,4)
      Call Gauss(nOrder+1,nOrder+1,B,A,FVal)
      XStart=Zero
      XLow = -One
      XHi  =  Two+Half
*
 500  Continue
c      Write(6,*) "After The call to Gauss", (A(I), I=0, 5)
c      Write(6,*) "Before the call  to Find_Min"
      Call Find_Min(nOrder,XStart,A,XMin,RC,XLow,XHi,ENew,Thr)
c      Write (*,*) 'ENew and Xmin=',ENew, Xmin
      If (.Not.RC) Then
*        Write (*,*) 'RC=.False.'
         Return
      End If
      If (Test2.lt.Zero) Then
         dqHdq=ENew-Energy(nIter)
      Else
         dqHdq=ENew
      End If
      XMin=XMin+Half
      Write (*,*) 'Minimum found at -->',XMin,'<--'
*                                                                      *
************************************************************************
*                                                                      *
*---- Update vectors in accordance with the result
*
*.... Compute the displacement iOld -> iNew
*
      Call DCopy(nInter,          q(1,iNew  ),1,dq(1,iOld),1)
      Call DaXpY(nInter,-One,     q(1,iOld  ),1,dq(1,iOld),1)
      dqdq=DDot(nInter,          dq(1,iOld),1,dq(1,iOld),1)
      Call DScal(nInter,XMin,    dq(1,iOld),1)
C      dqdq=DDot(nInter,          dq(1,iOld),1,dq(1,iOld),1)
*
*.... Compute the new q
*
      Call DCopy(nInter,          q(1,iOld  ),1, q(1,iNew  ),1)
      Call DaXpY(nInter,One,     dq(1,iOld),1,   q(1,iNew  ),1)
c      Print*, "The updated Geo.",  (q(I,1), I = 1, nInter)
c      Print*, "The line step", (dq(I,2), I = 1, nInter) 
*
*.... Update the gradient
*
      Call DScal(nInter,XMin,     g(1,iNew  ),1)
c      Print*, 'The gradients -1', (g(I, 1), I=1, nInter)
C      Call DaXpY(nInter,One-XMin, g(1,iOld  ),1, g(1,iNew  ),1)
      Call DaXpY(nInter, -XMin, g(1,iOld  ),1, g(1,iNew  ),1)
c      Print*, 'The gradients -2', (g(I, 1), I=1, nInter)
C      gdq=-DDot(nInter,           g(1,iNew  ),1,dq(1,iOld),1)
      gdq=DDot(nInter,           g(1,iNew  ),1,dq(1,iOld),1)

c      Write(6,*) 
c      Write(6,*) "scale factor", gdq/dqdq
      Call DaXpY(nInter,-gdq/dqdq,dq(1,iOld),1, g(1,iNew  ),1)
      Call DScal(nInter, -One, g(1, INew), 1)
c      Print*, 'The gradients', (g(I, 1), I=1, nInter)
*
*.... Update the displacement nIter-1 -> nIter
*
C      Call DCopy(nInter,          q(1,iNew  ),1,dq(1,iOld),1)
C      Call DaXpY(nInter,-One,     q(1,iOld),  1,dq(1,iOld),1)
*                                                                      *
************************************************************************
*                                                                      *
      Return 
      End
