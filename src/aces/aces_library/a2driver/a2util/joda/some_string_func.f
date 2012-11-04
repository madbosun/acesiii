      CHARACTER*4 FUNCTION some_string_func(NAME,Iord)
      CHARACTER*4 NAME, ItoA
      If (Name(2:2) .eq. 'N') then
         some_string_func(1:1) = Name(1:1)
         some_string_func(2: ) = ItoA(Iord, 0)
         some_string_func(linblnk(some_string_func)+1:) = Name(3:3)
      Else
         some_string_func = Name
      EndIf
      RETURN
      END
