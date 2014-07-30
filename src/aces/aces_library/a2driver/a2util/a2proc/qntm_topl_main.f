













































































































































































































































































































































































































































































































































      subroutine Qntm_Topl_Main(Icore, Icrsiz)
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
      logical bExist, Iwrite_H, Iwrite_T
      Character*4 Comp_pgrp, Full_pgrp
      Dimension Nocc(16), Atommass(Mxatms), Iatmchrg(Mxatms),  
     &          Fucoord(3,Mxatms), Coord(3,Mxatms),
     &          Norbits_fullG(Mxatms), NOrbits_compG(Mxatms),
     &          Nbsfns_4irrep(8)
      Dimension Icore(Icrsiz)
C
      Data Ione, Ieight, Iunit /1, 8, 10/
     
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
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Fucoord)
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
C
      print*, "Entering a2_prep4_mol_read"
      Print*, (Norbits_compG(i), i=1, Icmp_unq_atoms)
      Print*, (Norbits_fullG(i), i=1, Iful_unq_atoms)
      Call a2_prep4_mol_read(Icore, Maxcor, Natoms, Nbfns, 
     &                       Naobfns, Icmp_unq_atoms, 
     &                       Iful_unq_atoms, Coord, Iatmchrg, 
     &                       Norbits_fullG, Norbits_compG, 
     &                       Mxshel, Mxangmom, Itot_prim, Iunqshl, 
     &                       Ntotshl, Max_prim_4atom, 
     &                       Max_prim_4shell, Inext)

      Do i=1, 3
      Write(6, "(3F10.5)"), (coord(i,j), J=1,nreal_atoms)
      Enddo

      Write(6,*)
      Print*, "Offsets for various arrays initilized in Pre4_Den_plo."
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

C
      Icordens_a = Inext
      Icordens_b = Icordens_a + Nbfns*Nbfns*Iintfp
      Inatorbs_a = Icordens_b + Nbfns*Nbfns*Iintfp
      Inatorbs_b = Inatorbs_a + Naobfns*Nbfns*Iintfp
      Ieigvecs_a = Inatorbs_b + Naobfns*Nbfns*Iintfp
      Ieigvecs_b = Ieigvecs_a + Nbfns*Nbfns*Iintfp
      Ieigvals_a = Ieigvecs_b + Nbfns*Nbfns*Iintfp
      Ieigvals_b = Ieigvals_a + Nbfns*Iintfp
      Iscftmp1   = Ieigvals_b + Nbfns*Iintfp
      Iscftmp2   = Iscftmp1   + Nbfns*Nbfns*Iintfp 
      Iocca      = Iscftmp2   + Nbfns*Nbfns*Iintfp 
      Ioccb      = Iocca      + Nbfns*Iintfp
      Inext      = Ioccb      + Nbfns*Iintfp
      Write(6,*)
      print*, "Entering A2make_wfn, The geometry:"
    
      Do i=1, 3
      Write(6, "(3F10.5)") (coord(i,j), J=1,nreal_atoms)
      Enddo
      Call A2make_wfn(Icore(Icordens_a), Icore(Icordens_b),
     &                Icore(Inatorbs_a), Icore(Inatorbs_b),
     &                Icore(Ieigvecs_a), Icore(Ieigvecs_b),
     &                ICore(Ieigvals_a), Icore(Ieigvals_b), 
     &                Icore(Iscftmp1), Icore(Iscftmp2), 
     &                ICore(Iocca), Icore(Ioccb), Nocca, Noccb,
     &                Maxocca, Maxoccb, Naobfns, Nbfns, Iuhf,
     &                Tenergy)
C
       Lnp1             = Max_prim_4shell
       Icnterbf         = Inext 
       Inatorbs_reord_a = Icnterbf         + Nbfns
       Inatorbs_reord_b = Inatorbs_reord_a + Naobfns*Nbfns*Iintfp 
       Icoef_a          = Inatorbs_reord_b + Nbfns*Naobfns*Iintfp
       Icoef_b          = Icoef_a          + Nbfns*Itot_prim*Iintfp
       Itype            = Icoef_b          + Nbfns*Itot_prim*Iintfp
       ICent            = Itype            + Itot_prim*Iintfp
       IOrdr            = ICent            + Itot_prim*Iintfp
       IScr             = IOrdr            + Itot_prim
       INext            = Iscr             + Itot_prim*Iintfp
c
C
       Call Reord_mos(0, Nreal_atoms, Naobfns, Nbfns, Max_prim_4atom, 
     &                Icore(Icnterbf), ICore(Inatorbs_a),  
     &                ICore(Inatorbs_reord_a), Icore(Nangmom),
     &                Icore(Nmomao), Icore(Naoatm),.True.)
      Write(6,*)
      Write(6,*) "Alpha orbital for Quantum Topology analysis"
      Call output(Icore(Inatorbs_a), 1, nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
C
      CALL Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &            Icore(Ipcoeff),Itot_prim,Icore(Inatorbs_a),
     &            Naobfns,0.0D+00,Icore(Icoef_a),Itot_prim)

      Write(6,*) 
      Write(6,*) "Alpha orbital for Quantum Topology analysis"
      Call output(Icore(Icoef_a), 1, Itot_prim, 1, Nbfns, Itot_prim,
     &            Nbfns, 1)
      If (Iuhf .ne. 0) Then
         Call Reord_mos(0, Nreal_atoms, Naobfns, Nbfns, 
     &                  Max_prim_4atom, Icore(Icnterbf), 
     &                  ICore(Inatorbs_b), ICore(Inatorbs_reord_b), 
     &                  Icore(Nangmom), Icore(Nmomao), Icore(Naoatm),
     &                  .True.)
C
         Call Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &               Icore(Ipcoeff),Itot_prim,
     &               Icore(Inatorbs_b), Naobfns,0.0D+00,
     &               Icore(Icoef_b),Itot_prim)

      Write(6,*) 
      Write(6,*) "Beta orbital for Quantum Topology analysis"
      Call output(Icore(Icoef_b), 1, Itot_prim, 1, Nbfns, Itot_prim,
     &            Nbfns, 1)

      End if

      Call Nodummy(Natoms, Max_prim_4atom, Iatmchrg, Fucoord,
     &            Icore(Nangmom), Icore(Nmomfct), Icore(Inuctr))
C
      if (Comp_pgrp .Eq. "C1") Iful_unq_atoms = Icmp_unq_atoms
C
      Write(6,*)              
      Write(6,*) "The number of contracted function per shell"
      Write(6,*) (ICORE(ICONFUNTSHL - 1 + I), I=1, Ntotshl)
      Write(6,*) "The shell angular momentum"
      Write(6,*) (ICORE(IANGMOMTSHL - 1 + I), I=1, Ntotshl)
      Write(6,*) "The number of primitive functions per shell"
      Write(6,*) (ICORE(IPRMFUNTSHL - 1 + I), I=1, Ntotshl)
      Open(unit=Iunit, File="WFANAL", Form="Formatted", 
     &     Status="New")

      If (Iuhf .ne. 0) Then
         Iwrite_T = .False.
         Iwrite_H = .True.
      Else
         Iwrite_T = .True.
         Iwrite_H = .True.
      Endif 
C
      Call A2write_wfn(Icore(Icoef_a), Icore(Nangmom), Icore(Ialpha),
     &                 Icore(Iprmfuntshl), Icore(Iangmomtshl), 
     &                 Icore(Inprimshl), Icore(Itype), Icore(Icent),
     &                 ICore(Ieigvals_a), Icore(Iocca), Icore(Iordr),
     &                 ICore(Iscr),Coord, Maxocca, Mxshel, Itot_prim, 
     &                 Nbfns, Nreal_atoms,
     &                 Iatmchrg, Iunit, Iuhf, Iwrite_H, Iwrite_T, 
     &                 Tenergy)
C
      If (Iuhf .Eq. 0) Then
         Close (Iunit) 
      Else 
         Open(unit=Iunit, File="WFANAL", Form="Formatted",
     &        Status="Old")
         Ierr = 0
         Do While (Ierr .eq. 0)
            Read(Iunit, Iostat=Ierr) Temp_buff
         Enddo
      Endif
      If (Iuhf .ne. 0)  Call A2write_wfn(Icore(Icoef_b), 
     &                                   Icore(Nangmom), 
     &                                   Icore(Ialpha) , 
     &                                   Icore(Iprmfuntshl)   , 
     &                                   Icore(Iangmomtshl),
     &                                   Icore(Inprimshl),
     &                                   Icore(Itype), Icore(Icent),
     &                                   ICore(Ieigvals_b),
     &                                   Icore(Ioccb), Icore(Iordr),
     &                                   Icore(Iscr), Coord, Maxoccb, 
     &                                   Mxshel, Itot_prim, Nbfns, 
     &                                   Nreal_atoms, Iatmchrg,
     &                                   Iunit, Iuhf, .False., 
     &                                   .True., Tenergy)
C
      Return
      End

