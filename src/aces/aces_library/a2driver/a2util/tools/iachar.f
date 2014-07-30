
c This routine supplies the integer function 'iachar' if the Fortran
c compiler does not contain it as an intrinsic.

c Recall: iachar returns the decimal number of an ASCII character while
c         ichar returns the internal decimal of a character representation.


c This routine does NOTHING. Other files will include it so rigid
c Fortran compilers have something to do when a cpp define removes all
c of the compilable code in the .f file.

      subroutine to_combat_compiler_stupidity()
      return
      end





