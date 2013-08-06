#include "flags.h"
      subroutine A3_symadapt_main(Work, Icrsiz)
c
      implicit double precision (a-h,o-z)
c
#include "machsp.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "baslims.par"
C
      character*32 szFile
      logical bExist, Spherical
      Character*4 Comp_pgrp, Full_pgrp
      Dimension Atommass(Mxatms), Iatmchrg(Mxatms),  
     &          Fucoord(3,Mxatms), Coord(3,Mxatms),
     &          Norbits_fullG(Mxatms), NOrbits_compG(Mxatms),
     &          Nbfns_4irrep(8)
C
      Dimension Work(Icrsiz/iintfp)
C
      Data Ione, Ieight, Iunit /1, 8, 10/
C
      Iuhf = 1
      If (iflags(h_IFLAGS_ref).eq.0) iuhf = 0 
C
      Maxcor   = Icrsiz
      Length   = 0
c
c Read the JOBARC file for basic data of the molecule. 
c
      Call Getrec(20, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
      Call Getrec(20, 'JOBARC', 'NATOMS  ', Ione, Natoms)
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Fucoord)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
      Call Getrec(20, 'JOBARC', 'COMPNIRR', Ione, Nirrep)
      Call Getrec(20, 'JOBARC', 'NBASTOT ', Ione, Nbfns)
      Call Getrec(20, 'JOBARC', 'NAOBASFN', Ione, Naobfns)
      Call Getrec(20, 'JOBARC', 'NUMBASIR', Nirrep, Nbfns_4irrep)
C
      I = 1
      Do Iatm = 1, Natoms
         If (Iatmchrg(Iatm) .NE. 0) Then
            Coord(1, I) = Fucoord(1, Iatm)
            Coord(2, I) = Fucoord(2, Iatm)
            Coord(3, I) = Fucoord(3, Iatm)
            I = I + 1
         Endif
      Enddo
C
      Call Getcrec(20, 'JOBARC', 'COMPPTGP', 4, Comp_pgrp)
      Call Getcrec(20, 'JOBARC', 'FULLPTGP', 4, Full_pgrp)
      
      Spherical = .False. 
      If (iflags(h_IFLAGS_spherical).eq.1) Spherical = .True.

#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "Variable at symadapt main"
      Write(*, '(a)'), "Nreal_atoms, Natoms, Nirrep, Nbfns,"
      Write(*, '(a,5i4)')"Naobfns: ", Nreal_atoms, 
     &                    Natoms, Nirrep, Nbfns,Naobfns
      Write(6,*) "Spherical :", Spherical 
      Write(6,"(a,2A4)")"Comp_pgrp, Full_pgrp: ", Comp_pgrp, Full_pgrp
      Write(6,*)
#endif
C
      IScfvec_a = I0
      IScfvec_b = IScfvec_a  + Naobfns*Naobfns
      ITmp1     = IScfvec_b  + Naobfns*Naobfns
      ITmp2     = ITmp1      + Naobfns*Naobfns
      IScfevl_a = ITmp2      + Naobfns*Naobfns
      IScfevl_b = IScfevl_a  + Naobfns
      Ioed2a_sc = IScfevl_b  + Naobfns
      Ioed2a_or = Ioed2a_sc  + Naobfns
      Inext     = Ioed2a_or  + Naobfns

      Imemleft = (Icrsiz - Inext)
C
      Call a3_symadapt_scfvecs(Work(IScfvec_a), Work(IScfvec_b), 
     &                         Work(IScfevl_a), Work(IScfevl_b),
     &                         Work(Itmp1), Work(Itmp2),
     &                         Work(Ioed2a_sc), Work(Ioed2a_or),
     &                         Nbfns,
     &                         Naobfns, Nbfns_4irrep, Nirrep, Iuhf,
     &                         Spherical, Work(Inext), Imemleft)

      Return
      End

