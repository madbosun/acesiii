
c This routine supplies the character*1 function 'achar' if the Fortran
c compiler does not contain it as an intrinsic. This was first seen on
c an OSF/DEC compiler.

c Recall: achar returns the ASCII character of a decimal number while
c         char returns the internal character representation of a decimal.


c This routine does NOTHING. Other files will include it so rigid
c Fortran compilers have something to do when a cpp define removes all
c of the compilable code in the .f file.

      subroutine to_combat_compiler_stupidity()
      return
      end





