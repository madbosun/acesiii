#include "flags.h"
      subroutine Wrt_Prim_Main(Icore, Icrsiz)
c
      implicit double precision (a-h,o-z)
      parameter (mxcoef=30)
c
#include "machsp.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "baslims.par"
#include "offsets_2getden.com"
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
      Call Getrec(20, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
      Call Getrec(20, 'JOBARC', 'NATOMS  ', Ione, Natoms)
      Call Getrec(20, 'JOBARC', 'FULLNORB', Ione, Iful_unq_atoms)
      Call Getrec(20, 'JOBARC', 'COMPNORB', Ione, Icmp_unq_atoms)
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Coord)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
      Call Getrec(20, 'JOBARC', 'COMPNIRR', Ione, Nirrep)
      Call Getrec(20, 'JOBARC', 'OCCUPYA ', Nirrep, Nocc(1))
      Call Getrec(20, 'JOBARC', 'NBASTOT ', Ione, Nbfns)
      Call Getrec(20, 'JOBARC', 'NAOBASFN', Ione, Naobfns)
      Call Getrec(20, 'JOBARC', 'NUMBASIR', Nirrep, Nbsfns_4irrep)
      Call Getrec(20, 'JOBARC', 'FULLPOPV', Iful_unq_atoms, 
     &            Norbits_fullG)
      Call Getrec(20, 'JOBARC', 'COMPPOPV', Icmp_unq_atoms, 
     &            Norbits_compG)
C
      nOCCa = 0
      Do Irrep = 1, Nirrep
         Nocca = Nocca + Nocc(Irrep)
      End Do
      If (Iuhf .EQ. 1) Then
         Call Getrec(20, 'JOBARC', 'OCCUPYB', Nirrep, Nocc(I+8)) 
         Noccb = 0
         Do Irrep = 1, Nirrep
            Noccb = Noccb + Nocc(8+I) 
         End Do
      Else
         Noccb = Nocca
      Endif 
C
      Call Getcrec(20, 'JOBARC', 'COMPPTGP', 4, Comp_pgrp)
      Call Getcrec(20, 'JOBARC', 'FULLPTGP', 4, Full_pgrp)
C
#ifdef _DEBUG_LVLM2
      print*, "Entering A2_PREP4_MOL_READ.."
      Write(6,"(4I4)") (Norbits_compG(i), i=1, Icmp_unq_atoms)
      Write(6,*)
      Write(6,"(4I4)"), (Norbits_fullG(i), i=1, Iful_unq_atoms)
#endif 
      Call A2_PREP4_MOL_READ(Icore, Maxcor, Natoms, Nbfns, 
     &                       Naobfns, Icmp_unq_atoms, 
     &                       Iful_unq_atoms, Coord, Iatmchrg,
     &                       Norbits_fullG, Norbits_compG, 
     &                       Mxshel, Mxangmom, Itot_prim, Iunqshl, 
     &                       Ntotshl, Max_prim_4atom, 
     &                       Max_prim_4shell, Inext)

#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "Offsets for various arrays initilized in Wrt_Prim_Main."
      Write(*, '(31I5)'), IMEMBC,INUC,NUFCT,NUMOM,NFCT,NANGMOM,
     &                   IATMNAM,NAOUATM,NAOATM,IANGMOMSHL,
     &                   ICONFUNSHL,IPRMFUNSHL,INPRIMSHL,
     &                   IANGMOMTSHL,ICONFUNTSHL,IOFFSETPRM,
     &                   IOFFSETCON,IOFFSETSHL,IMAP_SHL2CNT,
     &                   IPRMFUNTSHL,NMOMAO,IALPHA,ICOEFFA,
     &                   ICOEFFB,IDENSA,IDENSB,IPCOEFFA,
     &                   IPCOEFFB, IPCOEFF,ISCR1,ISCR2
C
      Write(6,*) 
      Print*, "The orbital exponents in main", Itot_prim, ialpha
      call output(icore(ialpha), 1, Itot_prim, 1, 1, 1, Itot_prim,
     &            1,1,1)  
      Print*, "The Contraction Coefs."
      call output(icore(ipcoeff), 1, Itot_prim, 1, Naobfns, Itot_prim,
     &            Naobfns, 1)
#endif
       Itype            = Inext
       Icnt_xyz         = Itype            + Itot_prim
       Inext            = Icnt_xyz         + 3*Itot_prim*Iintfp
       Ileft            = Maxcor - Inext

       If (INext .ge. Maxcor) Call Insmem("Write_primt_main", Inext, 
     &                                     Maxcor)
C
      Call Write_prims(Natoms, Itot_prim, Nbfns, Naobfns, Coord, 
     &                 Icore(Ialpha), Icore(Ipcoeff), Ntotshl,
     &                 Icore(Inprimshl), Icore(Iangmomshl), 
     &                 Icore(Iconfunshl), Icore(Iprmfunshl),
     &                 Icore(Icnt_xyz), Icore(Itype), Iatmchrg)
C
      Return
      End

