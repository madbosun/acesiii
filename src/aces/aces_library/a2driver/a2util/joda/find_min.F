      Subroutine Find_Min(nOrder,XStart,A,XMin,RC,XLow,XHi,ENew,
     &                    THR)
      Implicit Double Precision (a-h,o-z)

      DIMENSION A(0:nOrder)
      Logical RC
      DATA ZERO /0.0D0/, ONE /1.0D0/
*
      XValue=XStart
      RC=.True.
c      Write(6,*) nOrder,XStart,XMin,RC,XLow,XHi,ENew,
c     &                    THR
      MaxIter=100
      Do i = 1, MaxIter
         X=XValue
         fnc  = Zero
         XX   = One 
         Do j=0,nOrder
            fnc=fnc+A(j)*XX
            XX = XX*X
         End Do
         dfnc = Zero
         XX   = One 
         Do j=1,nOrder
            tmp = j
            dfnc=dfnc+A(j)*tmp*XX
            XX = XX*X
         End Do
         ddfnc = Zero
         XX    = One 
         Do j=2,nOrder
            tmp = j*j-j
            ddfnc=ddfnc+A(j)*tmp*XX
            XX = XX*X
         End Do
         XInc=dfnc/ddfnc
         XValue=XValue-XInc
c
c         Write (6,*) 'Fnc,dFnc,ddFnc,xinc=',Fnc,dFnc,ddFnc,xinc,
c     &                xvalue
         If (Abs(XInc).lt.Thr) Then
            ENew=fnc
            XMin=XValue
C            Call QExit('Find_Min')
            Return
         End If
         If (XValue.gt.XHi) XValue=XHi
         If (XValue.lt.XLow) XValue=XLow
      End Do
      Write (*,*) '-- Too many iterations in Find_Min'
      RC=.False.
*
      Return
      End

