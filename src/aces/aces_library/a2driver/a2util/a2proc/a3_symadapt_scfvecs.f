      Subroutine  a3_symadapt_scfvecs(Scfvec_a, Scfvec_b, Scfevl_a,
     &                                Scfeval_b, Tmp1, Tmp2, 
     &                                Oed2AScal, Ioed2Aord, 
     &                                Nbfns,
     &                                Naobfns, Nbfirr, Nirrep, Iuhf,
     &                                Spherical, Work, Imemleft)

      Implicit Double Precision (A-H, O-Z)
      Logical Spherical
      Parameter (Tol = 1.0D-09)



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




      Dimension Scfvec_a(Naobfns*Naobfns), Scfvec_b(Naobfns*Naobfns), 
     &          Scfevl_a(Nbfns), Scfeval_b(Nbfns), 
     &          Tmp1(Naobfns*Naobfns), Nbfirr(8),
     &          Tmp2_a(Naobfns*Naobfns), Tmp2_b(Naobfns*Naobfns), 
     &          Oed2AScale(Naobfns), Ioed2Aorder(Naobfns)

      Dimension Nocc(8,2)
     
      Call Getrec(20, "JOBARC", "SCFEVCA0", Nbfns*Nbfns*Iintfp,
     &            Tmp2_a)
      Call Filter(Tmp2_a, Nbfns*Nbfns, Tol)
      Call Getrec(20, "JOBARC", "SCFEVLA0", Nbfns*Iintfp, Scfevl_a)

      If (Iuhf .EQ. 1) then
         Call Getrec(20, "JOBARC", "SCFEVCB0", Nbfns*Nbfns*Iintfp,
     &               Tmp2_b)
      Call Filter(Tmp2_b, Nbfns*Nbfns, Tol)
         Call Getrec(20, "JOBARC", "SCFEVLB0", Nbfns*Iintfp, Scfevl_b)

      Endif 

      Write(6,"(a)") "The SCF eigenvectors from ACES III run"
      Call output(Tmp2_a, 1, Nbfns, 1, Nbfns, Nbfns, Nbfns, 1)
      Write(6,"(a)") "The Eigenvalues"
      Write(6,"(6(1x,F10.5))") (Scfevl_a(i), i=1, Nbfns)
      if (Iuhf .Gt. 0) Call output(Tmp2_b, 1, Nbfns, 1, Nbfns, 
     &                             Nbfns, Nbfns, 1)

C Get the ERD to ACES scaling and ordering vectors.

      Call Getrec(20, "JOBARC", "ERD2A2CS", Nbfns*Iintfp,
     &            Oed2AScale)
      Call Getrec(20, "JOBARC", "ERDORDER", Nbfns, Ioed2Aorder)
C
      Call Do_oed_to_vmol(Nbfns, Ioed2Aorder, Oed2AScale, Tmp2_a,
     &                    Scfvec_a)
      If (Iuhf .EQ. 1) then
         Call Do_oed_to_vmol(Nbfns, Ioed2Aorder, Oed2AScale, Tmp2_b,
     &                      Scfvec_b)
      Endif
C 
C Generate the occupation numbers for each irrep based on eigen
C values and the number of basis functions per irrep.

      Call Occupy(Nirrep, Nbfirr, Nbfns, Scfevl_a, Work, Nocc(1,1),
     &            1)
      If (Iuhf .EQ. 1) Then
         Call Occupy(Nirrep, Nbfirr, Nbfns, Scfevl_b, Work, Nocc(1,2),
     &               2)
      Else
         Call Icopy(8, Nocc(1, 1), 1, Nocc(1,2), 1)
      Endif
C
      Write(6,"(a)") "The occupation numbers"
      Write(6,"(8(1x,I3))") (Nocc(i, 1), i=1, Nirrep)
      Write(6,"(8(1x,I3))") (Nocc(i, 2), i=1, Nirrep)
C
C First convert from Spherical to Cartesian (if the calculation is
C in Cartesian this should do nothing).

      Call Getrec(20, "JOBARC", "CMP2ZMAT", Nbfns*Naobfns*Iintfp,
     &            Tmp1)

      Call Xgemm("N", "N", Naobfns, Nbfns, Nbfns, 1.0D0, Tmp1,
     &            Naobfns, Scfvec_a, Nbfns, 0.0D0, Tmp2, Naobfns)

      Call Dcopy(Naobfns*Nbfns, Tmp2, 1, Scfvec_a, 1)

      If (Iuhf .EQ. 1) Then
         Call Xgemm("N", "N", Naobfns, Nbfns, Nbfns, 1.0D0, Tmp1, 
     &               Naobfns, Scfvec_b, Nbfns, 0.0D0, Tmp2, Naobfns)

         Call Dcopy(Naobfns*Nbfns, Tmp2, 1, Scfvec_b, 1)

      Endif

      Write(6,*) "The SCF eigenvectors from ACES III (Cartesian basis)"
      Call output(Scfvec_a, 1, Naobfns, 1, Nbfns, Naobfns, Nbfns, 1)
      if (Iuhf .Gt. 0) Call output(Scfvec_b, 1, Naobfns, 1, Nbfns, 
     &                             Naobfns, Nbfns, 1)

      Call Get_irreps(Scfvec_a, Scfevl_a, Work, Imemleft*Iintfp, 
     &                Nbfns, Naobfns, 1, Nocc, Iuhf)
      If (Iuhf .EQ. 1) Call Get_irreps(Scfvec_b, Scfevl_b, Work, 
     &                                 Imemleft*Iintfp, Nbfns, 
     &                                 Naobfns, 2, Nocc, Iuhf)
C
      Return
      End

