      Subroutine PitPtG (RealPG, PitS, PitOrd)
C
C     Take a real Schoenflies symbol & convert it into something Pitzer
C     can understand.  Pretty inelegant for now.  I hope we won't need
C     to use Pitzer much longer.
C
      Character*(*) RealPG, PitS
      Integer PitOrd
      Integer AtoI
C
      Do 100 i = 1, Len(PitS)
         PitS(i:i) = ' '
 100  Continue
      If (RealPG(1:linblnk(RealPG)) .eq. 'C1') then
         PitS = 'C1'
         PitOrd = 0
         Return
      EndIf
      If (RealPG(1:linblnk(RealPG)) .eq. 'CS' .OR.
     $   RealPG(1:linblnk(RealPG)) .eq. 'Cs' .OR.
     $   RealPG(1:linblnk(RealPG)) .eq. 'C S' .OR.
     $   RealPG(1:linblnk(RealPG)) .eq. 'C s') then
         PitS = 'CS'
         PitOrd = 0
         Return
      EndIf
      If (RealPG(1:1) .eq. 'C' .OR. RealPG(1:1) .eq. 'D') then
         PitS(1:1) = RealPG(1:1)
         PitS(2:2) = 'N'
         PitOrd = AtoI(RealPG(2:))
         i = linblnk(RealPG)
         If (RealPG(i:i) .eq. 'v' .OR. RealPG(i:i) .eq. 'V') then
            PitS(3:3) = 'V'
         ElseIf (RealPG(i:i) .eq. 'h' .OR. RealPG(i:i) .eq. 'H') then
            PitS(3:3) = 'H'
         EndIf
      Else
         PitS = 'XXX'
      EndIf
      Return
      End
