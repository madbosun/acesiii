













































































































































































































































































































































































































































































































































      subroutine Den_Plots_Main(Icore, Icrsiz)
c
      implicit double precision (a-h,o-z)
      parameter (mxcoef=30)
c


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end










c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)


      COMMON /OFF_4INTGRT/IMEMBC,INUC,NUFCT,NUMOM,NFCT,NANGMOM,
     &                    IATMNAM,NAOUATM,NAOATM,IANGMOMSHL,
     &                    ICONFUNSHL,IPRMFUNSHL,INPRIMSHL,
     &                    IANGMOMTSHL,ICONFUNTSHL,IOFFSETPRM,
     &                    IOFFSETCON,IOFFSETSHL,IMAP_SHL2CNT,
     &                    IPRMFUNTSHL,NMOMAO,IALPHA,ICOEFFA,
     &                    ICOEFFB,IDENSA,IDENSB,IPCOEFFA,
     &                    IPCOEFFB, IPCOEFF,ISCR1,ISCR2,IREORD 




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
      If (iflags(11).eq.0) iuhf = 0 
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
      Call A2_PREP4_MOL_READ(Icore, Maxcor, Natoms, Nbfns, 
     &                       Naobfns, Icmp_unq_atoms, 
     &                       Iful_unq_atoms, Coord, Iatmchrg,
     &                       Norbits_fullG, Norbits_compG, 
     &                       Mxshel, Mxangmom, Itot_prim, Iunqshl, 
     &                       Ntotshl, Max_prim_4atom, 
     &                       Max_prim_4shell, Inext)

C
       Lnp1             = Itot_prim
       Icnterbf         = Inext 
       Iscfcoef_a       = Icnterbf         + Nbfns
       Iscfcoef_b       = Iscfcoef_a       + Naobfns*Nbfns*Iintfp
       Iscfcoef_reord_a = Iscfcoef_b       + Naobfns*Nbfns*Iintfp 
       iSCfcoef_reord_b = Iscfcoef_reord_a + Nbfns*Naobfns*Iintfp
       Icoef_a          = Iscfcoef_reord_b + Nbfns*Naobfns*Iintfp
       Icoef_b          = Icoef_a          + Nbfns*Itot_prim*Iintfp
       Itmp1            = Icoef_b          + Nbfns*Itot_prim*Iintfp
       Itmp2            = Itmp1            + Lnp1*Lnp1*Iintfp
       Itmp3            = Itmp2            + Lnp1*Lnp1*Iintfp
       Iprductint       = Itmp3            + Lnp1*Lnp1*Iintfp
       Inuctr           = Iprductint       + Naobfns*Naobfns*Iintfp
       Itype            = Inuctr           + Natoms
       Icnt_xyz         = Itype            + Itot_prim
       Inext            = Icnt_xyz         + 3*Itot_prim*Iintfp
       Ileft            = Maxcor - Inext

       If (INext .ge. Maxcor) Call Insmem("Den_plot_main", Inext, 
     &                                     Maxcor)
C
      print*, "Entering Reord_MOs:"
      Print*, "Input variables", Inext, Max_prim_4atom, Itot_prim, 
     &         icnterbf, Iscfcoef_reord_a, Iscfcoef_reord_b,
     &         Icoef_a, Icoef_b, Nangmom, Nmomao, Naoatm,ipcoeff
      write(*,*)
C
       Call Reord_mos(0, Nreal_atoms, Naobfns, Nbfns, Max_prim_4atom, 
     &                Icore(Icnterbf), ICore(Iscfcoef_a),  
     &                ICore(Iscfcoef_reord_a), Icore(Nangmom),
     &                Icore(Nmomao), Icore(Naoatm), .False.)
C
      CALL Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &            Icore(Ipcoeff),Itot_prim,Icore(Iscfcoef_reord_a),
     &            Naobfns,0.0D+00,Icore(Icoef_a),Itot_prim)


      If (Iuhf .ne. 0) Then
         Call Reord_mos(1, Nreal_atoms, Naobfns, Nbfns, 
     &                  Max_prim_4atom, Icore(Icnterbf), 
     &                  ICore(Iscfcoef_b), ICore(Iscfcoef_reord_b), 
     &                  Icore(Nangmom), Icore(Nmomao), Icore(Naoatm),
     &                  .False.)
C
         Call Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &               Icore(Ipcoeff),Itot_prim,
     &               Icore(Iscfcoef_reord_b), Naobfns,0.0D+00,
     &               Icore(Icoef_b),Itot_prim)
      End if

      print*, "Entering Nodummy:"
      write(*,*)
      CALL Nodummy(Natoms, Max_prim_4atom, Iatmchrg, Coord,
     &            Icore(Nangmom), Icore(Nmomfct), Icore(Inuctr))
C
      if (Comp_pgrp .Eq. "C1") Iful_unq_atoms = Icmp_unq_atoms
C
      Call Generate_12Dgrid(Icore(Inext), Ileft, Natoms, Itot_prim,
     &                      Nbfns, Naobfns, Lnp1, Norbits_compG,
     &                      Iful_unq_atoms, Coord, Icore(Ialpha),
     &                      Icore(Ipcoeff), Ntotshl,Icore(Inprimshl), 
     &                      Icore(Iangmomshl), Icore(Iconfunshl), 
     &                      Icore(Iprmfunshl), Icore(Nufct), 
     &                      Icore(Naouatm), Icore(Ioffsetprm),
     &                      Icore(Ioffsetcon), Icore(Itmp1),
     &                      Icore(Itmp2), Icore(Itmp3),
     &                      Icore(Iprductint), Icore(Icnt_xyz),
     &                      Icore(Ireord), Icore(Itype), 
     &                      Iatmchrg, Iuhf)
      write(6,*)
      Write(6,"(t3,50a)")"@-den_plot_main: The density distribution ",
     &                  "calculation is succesfully completed."
      Return
      End

