      Subroutine Post_opt_update
c
      Implicit Double Precision (A-H, O-Z)
c
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)




































































































































































































































































































































































































































































































































c     Maximum string length of basis set
      INTEGER BASLEN
      PARAMETER (BASLEN=80)







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



C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end


c
      Integer Genby(Mxatms) 
      Double Precision Orient(3, 3)
      Character*(baslen) BasNam(Mxatms)
      Character*4 FPGrp, BPGrp, PGrp 
      Character*5, Zsym(Mxatms)
      Dimension Scratch(Mxatms) 
      Logical Can_do_freq 
c
      Common /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
c
      Call Entry(BasNam, .True., Can_do_freq)
c
      Call Getrec(20, 'JOBARC', 'ZMATATMS', 1, NATOMS)
      Call Getrec(20, 'JOBARC', 'LINEAR  ', 1, ILINEAR)
      If (ILINEAR .EQ. 1) Then
          NX   = 3*NATOMS
          NXM6 = NX - 5
      Else
          NX   = 3*NATOMS
          NXM6 = NX -6
      Endif  
c
      Call Getrec(20, "JOBARC", "COORD   ", NX*IINTFP, Q)
      Call Getrec(20, 'JOBARC', 'CONCTVTY', NX, NCON)
      Call Getrec(-1, 'JOBARC', 'CORD_INT', NX*IINTFP, R)
      Call Bohr2angs(R, NX)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', NATOMS*IINTFP,
     &           ATMASS)
      Call Getrec(20, 'JOBARC', 'ORIENTMT', 9*IINTFP,
     &            ORIENT)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', NATOMS, IATNUM)
      Call Getrec(20, 'JOBARC', 'ICSQUASH', NX, ISQUASH)
      Call Getcrec(20, 'JOBARC', "PTGP    ", 4, FPGRP)
      Call Getcrec(20, 'JOBARC', "ABL_PTGP", 4, BPGRP)
      Call Getcrec(20, 'JOBARC', "CMP_PTGP", 4, PGRP)
      Call Getcrec(20, 'JOBARC', "ZSYM", 5*NATOMS, ZSYM)
           Print*, "The data read from JOBARC in post_opt_update"
           Print*, "The NATOMS:", NATOMS
           Print*, "Internal coords:",(R(I), I=1, NX)
           Print*, "The Cartesian coords:", (Q(I), I=1, NX) 
           print*, "The connectivities:", (NCON(I), I=1, NX)
           Print*, "The atomic charges:", (IATNUM(I), I=1, NX/3)
           Print*, "The point group", PGRP 
c
c The scratch and Nunique are not used. Genby, istat and Genby
c are internal 

      Call MkVMOL(Q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &            GenBy, Scratch, IStat, BasNam)
      
      Print*, "out from Mkvmol in post_opt_update"
      Return
      End


