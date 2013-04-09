#include "flags.h"
      subroutine Read_File_Main(Icore, Icrsiz)
c
      implicit double precision (a-h,o-z)
      parameter (mxcoef=30)
c
#include "machsp.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "baslims.par"
C
      character*32 szFile
      logical bExist
      Character*4 Comp_pgrp, Full_pgrp
      Dimension Nocc(16), Atommass(Mxatms), Iatmchrg(Mxatms),  
     &          Coord(3*Mxatms), Norbits_fullG(Mxatms), 
     &          NOrbits_compG(Mxatms),Nbsfns_4irrep(8)
      Dimension Icore(Icrsiz)
C
      Data Ione, Ieight /1, 8/
C
      Iuhf = 1
      If (iflags(h_IFLAGS_ref).eq.0) iuhf = 0 
C
      Maxcor   = Icrsiz
      Mxangmom = Nht
      Length   = 0
c
c Read the JOBARC file for basic data of the molecule. 
c
      
      Call Getrec(-1, 'JOBARC', 'NUMDROPA', Length, Ijunk)       
      If (Length .GT. 0) Then 
         Print*, "Frozen-core is not allowed density plots"
         Call Errex
      Endif
c 
      Call Getrec(-1, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
      Call Getrec(-1, 'JOBARC', 'NATOMS  ', Ione, Natoms)
      Call Getrec(-1, 'JOBARC', 'FULLNORB', Ione, Iful_unq_atoms)
      Call Getrec(-1, 'JOBARC', 'COMPNORB', Ione, Icmp_unq_atoms)
      Call Getrec(-1, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Coord)
      Call Getrec(-1, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
      Call Getrec(-1, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
      Call Getrec(-1, 'JOBARC', 'COMPNIRR', Ione, Nirrep)
      Call Getrec(-1, 'JOBARC', 'OCCUPYA ', Nirrep, Nocc(1))
      Call Getrec(-1, 'JOBARC', 'NBASTOT ', Ione, Nbfns)
      Call Getrec(-1, 'JOBARC', 'NAOBASFN', Ione, Naobfns)
      Call Getrec(-1, 'JOBARC', 'NUMBASIR', Nirrep, Nbsfns_4irrep)
      Call Getrec(-1, 'JOBARC', 'FULLPOPV', Iful_unq_atoms, 
     &            Norbits_fullG)
      Call Getrec(-1, 'JOBARC', 'COMPPOPV', Icmp_unq_atoms, 
     &            Norbits_compG)
C
       Nvibs  = (3*Nreal_atoms- 6) 
       Iloc   = 1
       Iscr1  = Iloc   + 9*Nvibs*Nvibs
       INext  = IScr1  + 9*Nvibs*Nvibs
       Ilef   = Maxcor - Inext

       If (INext .ge. Maxcor) Call Insmem("Write_primt_main", Inext, 
     &                                     Maxcor)
C
       Call Read_filE(Nreal_atoms, Nvibs, Icore(Iloc), ICore(IScr1))
C
      Return
      End

