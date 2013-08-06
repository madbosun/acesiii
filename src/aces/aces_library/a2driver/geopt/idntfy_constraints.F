      Subroutine Idntfy_constraints(Iredunco, Nredunco, TOTNOFBND, 
     &                              TOTNOFANG, TOTNOFDIH, Not_found)
C
      Implicit Double Precision (A-H, O-Z)
C
#include "mxatms.par" 
C
      Integer Constrains, Total_constrains, TOTNOFBND, TOTNOFANG,
     &        TOTNOFDIH, Constrained
      Logical Not_found, Sucess
C
      Dimension Iredunco(4, Nredunco), Constrains(4, Maxredunco),
     &          Constrained(Maxredunco)
C 
C Read the constrins from the *CONSTRAINS namelist from the
C ZMAT file.
C
#ifdef _DEBUG_LVL0
      write(6,"(a)") " Entering read_constrains"
#endif
      Call Read_ric_constrains(Constrains,
     &                         Num_bond_Constrains,
     &                         Num_angle_constrains,
     &                         Num_dihs_angle_constrains,
     &                         Not_found)

      If (Not_found) Return

#ifdef _DEBUG_LVL0
      Write(6,*)
      write(6,"(a,3(1x,I4))") " Constrains:bonds,angles,dihedrals ",
     &                        Num_bond_Constrains, 
     &                        Num_angle_Constrains,
     &                        Num_dihs_angle_Constrains
#endif
      Total_constrains = Num_bond_Constrains + Num_angle_Constrains +
     &                   Num_dihs_angle_Constrains

#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The constrains read from ZMAT"
      Write(6,*)
      Do i = 1, Total_constrains
         Write(6,111) (Constrains(j, i), j=1, 4)
      Enddo
 111  Format(5X, 4(I3, 1X))
#endif
C
C If there are no constrains to be found, return and do the full
C optimization. 
C

      Sucess = .False. 
      Ifound = 0
      Call Izero(Constrained, Nredunco)
C
      Do Jc = 1, Num_bond_Constrains
         Do Ic = 1, TOTNOFBND
            
            ICcon = Constrains(1, Jc)
            JCcon = Constrains(2, Jc)

            IRcon = Iredunco(1, Ic)
            JRcon = Iredunco(2, Ic)
#ifdef _DEBUG_LVLM2
            Write(6,"(4i5)") ICcon, JcCon, IRcon, JRcon
#endif
            If (ICcon .eq. IRcon .OR. ICcon .eq. JRcon) Then

                If (JCcon .Eq. IRcon .OR. JCcon .eq. JRcon) Then

                    Constrained(Ic) = 1
                    Ifound = Ifound + 1
                    Sucess = .TRUE.

                Endif

            Endif 
            IF (Sucess) GO TO 10

        Enddo
  10    Continue 
        Sucess = .False. 
      Enddo
            
      Ibegin = Num_bond_Constrains + 1
      Jbegin = TOTNOFBND + 1

      Do Jc = Ibegin, Ibegin + Num_angle_Constrains - 1
         Do Ic = Jbegin, Jbegin + TOTNOFANG - 1
            
            ICcon = Constrains(1, Jc) 
            JCcon = Constrains(2, Jc)
            KCcon = Constrains(3, Jc)

            IRcon = Iredunco(1, Ic)
            JRcon = Iredunco(2, Ic)
            KRcon = Iredunco(3, Ic)

#ifdef _DEBUG_LVLM2
            Write(6,"(6i5)") ICcon, JcCon, KCcon, IRcon, JRcon, KRcon
#endif

            If (ICCon .eq. IRcon .OR. ICCon .eq. JRcon .OR. 
     &          ICcon .EQ. KRcon) Then
                If (JCCon .eq. IRcon .OR. JCCon .eq. JRcon .OR.
     &              JCcon .EQ. KRcon) Then
                    If (KCCon .eq. IRcon .OR. KCCon .eq. JRcon .OR.
     &                  KCcon .EQ. KRcon) Then

                        Constrained(Ic) = 1
                        Ifound = Ifound + 1
                        Sucess = .TRUE.

                     Endif
                Endif
            Endif

            IF (Sucess) GO TO 20

         Enddo
  20     Continue 
         Sucess = .False. 
      Enddo
C
#ifdef _DEBUG_LVLM
      Write(6,*) 
      Write(6,*) "The CONSTRAINED array"
      Write(6,"(6(1x,I4))") (Constrained(i),i=1,Nredunco)
#endif
      Ibegin = Num_bond_constrains + Num_angle_constrains + 1
      Jbegin = TOTNOFBND + TOTNOFANG + 1

      Do Jc = Ibegin, Ibegin + Num_dihs_angle_Constrains - 1
         Do Ic = Jbegin, Jbegin + TOTNOFDIH - 1

            ICcon = Constrains(1, Jc)
            JCcon = Constrains(2, Jc)
            KCcon = Constrains(3, Jc)
            LCcon = Constrains(4, Jc)

            IRcon = Iredunco(1, Ic)
            JRcon = Iredunco(2, Ic)
            KRcon = Iredunco(3, Ic)
            LRcon = Iredunco(4, Ic)
#ifdef _DEBUG_LVLM2
            Write(6,"(8i5)") ICcon, JcCon, KCcon, LCcon, IRcon, JRcon, 
     &                 KRcon,LRcon
#endif

            If (ICCon .eq. Ircon .OR. ICcon .Eq. JRcon .OR.
     &          ICcon .eq. Krcon .OR. ICCon .EQ. LRcon) Then
                If (JCCon .eq. Ircon .OR. JCcon .Eq. JRcon .OR.
     &              JCcon .eq. Krcon .OR. JCCon .EQ. LRcon) Then
                    If (KCCon .eq. Ircon .OR. KCcon .Eq. JRcon .OR.
     &                  KCcon .eq. Krcon .OR. KCCon .EQ. LRcon) Then
                        If (LCCon .eq. Ircon .OR. LCcon .Eq. JRcon .OR.
     &                      LCcon .eq. Krcon .OR. LCCon .EQ. LRcon) Then

                            Constrained(Ic) = 1
                            Ifound = Ifound + 1
                            Sucess = .TRUE.

                        Endif 
                    Endif
                Endif
            Endif

            IF (Sucess) GO TO 30

         Enddo
  30     Continue
         Sucess = .False.
      Enddo
C
      If (Ifound .Lt. Total_constrains) Then
        Write(6,"(3x,2a)")"At least one constrain specified" ,
     &               " in the input does not correspond to a known"
        Write(6,"(3x,2a)") "RIC. The constrain(s) is(are) ignored",
     &               " and the optimization proceeds."
C
        Write(6,*)
        Write(6,"(3x,2a)") "If this has been a mistake, cancel the job",
     &                " and use the redundent coordinates"
        Write(6,"(3x,a,a)")"printed below to identify the correct ", 
     &                     "degs. of freedoms that can be " 
        Write(6,"(3x,a)")  "constrained."
C
        Write(6,*)
        Write(6,"(3x,a)") "The center lables of redundent coordinates"
        Do i = 1, Nredunco
           Write(6,112) (iredunco(j, i), j=1, 4)
        Enddo
 112    Format(5X, 4(I3, 1X))

      Endif
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) "The CONSTRAINED array"
      Write(6,"(a,I4)") "The total No. of constraints found",Ifound
      Write(6,"(6(1x,I4))") (Constrained(i),i=1,Nredunco)
#endif
C
      Call IPutrec(20, "JOBARC","CONSTRNS", Nredunco, Constrained)
C
      Return
      End
        
