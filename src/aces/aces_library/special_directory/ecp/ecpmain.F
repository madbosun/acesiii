#include "flags.h"
C
      Subroutine ecpmain() 
      
      Implicit Double Precision (A-H, O-Z)

#include "ecp.par"
#include "ecp.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "icore.com"

      Character*4 INP_Fname
      Character*4 ATMNAM(Max_centers), Namat(Max_centers)
    
      Logical Cartesian, Grads, Spherical
      Integer Angmom, Atomic_label
      Integer*8 Istart_Heap, IECP_int, IECP_Grad_int
  
      DIMENSION WORK(Maxmem_a3), Decp_ints(1)
    
      DIMENSION NFCT(max_centers), NAOATM(max_centers),
     &          NUFCT(max_centers), NAOUATM(max_centers),
     &          NAOUATM2(max_centers),
     &          iqmstr(max_centers),jcostr(max_centers,Maxjco),
     &          nucstr(max_centers,Maxang,Maxjco),
     &          nrcstr(max_centers,Maxang,Maxjco),
     &          Jstrt(max_shells),Jrs(max_shells),
     &          NSpc(max_centers), NCFpS(max_shells), 
     &          NPFpS(max_shells), Atomic_label(max_shells),
     &          Charge(max_centers),IGenby(max_centers),
     &          IXALPHA(max_shells),Angmom(Max_shells),
     &          NPOP(max_centers),IXPCOEF(max_shells),
     &          IREORDER(max_centers),
     &          COORD(3,max_shells),Cord(3,max_centers)
    
C
      If (IFLAGS(h_IFLAGS_ecp) .EQ. 1) Iecp = 1
      If (IFLAGS(h_IFLAGS_spherical) .EQ. 1) Cartesian = .True. 
      If (IFLAGS2(h_IFLAGS2_grad_calc) .EQ. 1) Grads = .True. 

      INP_Fname = "MOL"
      Call SIMPLE_INSPECT_MOLECP(INP_Fname, max_centers, max_shells,
     &                        ncenters, nshells, nspc, cartesian, 
     &                        ITFCT, LNP1, lnpo, nfct, nufct, 
     &                        nbasis, NAOBASIS, nCFpS, 
     &                        nPFpS, NAOATM, angmom, atomic_label, 
     &                        vnn, Maxang, Maxjco, Iecp,
     &                        NUcstr, Nrcstr, Iqmstr, Jcostr, 
     &                        Jstrt, Jrs, Atmnam, Charge, NAOUATM,
     &                        NAOUATM2)

C      Write(6,*) "The NJCO", NJCO
C---
c      Write(6,"(a,a,7(1x,I4))") "natoms,nshells,ITFCT,LNP1,lnpo",
c     &                          "nbasis,NAOBASIS:",ncenters,nshells,
c     &                                         ITFCT,LNP1,lnpo,
c     &                                        nbasis,NAOBASIS
c      Write(6,*) 
c      Write(6,"(a)") "The number of con. function per shell"
c      Write(6,"(6(1x,I4))") (nCFpS(i), i=1, nshells)
c      Write(6,"(a)") "The number of prim. function per shell"
c      Write(6,"(6(1x,I4))") (nPFpS(i), i=1, nshells)
c      Write(6,"(a)") "The Iqmstr"
c      Write(6,"(6(1x,I4))") (Iqmstr(i), i=1, ncenters)
c      Write(6,"(a)") "The jcostr"
c      Do i=1, Ncenters
c        Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
c      Enddo
c      Write(6,"(a)") "The Nucstr"
c      Do i=1, Ncenters
c        Write(6,"(4(1x,I4))") ((Nucstr(i,j,k),j=1,4),k=1,Maxjco)
c      Enddo
c      Write(6,"(a)") "The Nrcstr"
c      Do i=1, Ncenters
c        Write(6,"(4(1x,I4))") ((Nrcstr(i,j,k),j=1,4),k=1,Maxjco)
c      Enddo
c      Write(6,"(a)") "The jstrt"
c      Write(6,"(6(1x,I4))") (Jstrt(i), i=1, Nshells)
c      Write(6,"(a)") "The jrs"
c      Write(6,"(6(1x,I4))") (Jrs(i), i=1, Nshells)

      nalpha = 0
      npcoef = 0
      do i = 1, nshells
         nalpha = nalpha + npfps(i)
         npcoef = npcoef + npfps(i) * ncfps(i)
      enddo

      IALPHA   = 1
      IPCOEF   = IALPHA +  nalpha
      ISALPHA  = IPCOEF +  npcoef
      ISPCOEF  = ISALPHA + LNP1
      INEXT    = ISPCOEF + LNPO

      do i = 1, ncenters
         npop(i)     = 1
         ireorder(i) = i
         Namat(i) = ATMNAM(i)
      enddo
  
      If (Inext .GE. Maxmem_a3)  Then
          Write(6, "(1x,a)") "@-ecpmain: insufficent memory."
          Write(6, "(1x,a,1x,I5)") "Needed (doubles):", Inext
          Write(6, "(1x,a,1x,I5)") "Available (doubles):", Maxmem_a3
          Call Errex
      Endif
 
      Call READ_BASIS_INFOECP(INP_FNAME, Ncenters, Ncenters, NPOP,
     &                     IREORDER,
     &                     CARTESIAN, ITFCT, LNP1, LNPO,
     &                     NFCT, NBASIS, WORK(IALPHA), IXALPHA,
     &                     Work(IPCOEF), IXPCOEF, Work(ISALPHA), 
     &                     Work(ISPCOEF), MAX_CENTERS, ATMNAM, 
     &                     COORD, CORD, NAOATM)

c     Write(6,*)
c     Write(6,"(a)") "The exponents"
c     Write(6,"(6(1x,F12.8))") (Work(Ialpha+i), i=0, Nalpha-1)
c     Write(6,*)
c     Write(6,"(a)") "The contraction coefs"
c     Write(6,"(6(1x,F12.8))") (Work(Ipcoef+i), i=0, Npcoef-1)
c     Write(6,"(a)") "The jcostr"
c     Do i=1, Ncenters
c        Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
c     Enddo

c     Write(6,*)
c     Write(6,*) "Atom Coordiantes of each shell"
c     Do Iatm =1, Ncenters
c        Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
c     Enddo
   
      IECP_int      = 1
      IECP_grad_int = IECP_int +  (1/2)*(NAOBASIS+1)*NAOBASIS

      If (Grads) Then
          Mem_needed    = IECP_grad_int + (3/2)*(NAOBASIS+1)*NAOBASIS*
     &                                     Ncenters
      Else
          Mem_needed    = IECP_grad_int
      Endif

      Ierr = 0
      Call MEM_ALLOC(Decp_ints, Mem_needed, Ifltln, Istart_Heap, .True.,
     &               Ierr)
C
      If (Ierr .ne. 0) then
         Write(6, "(a,i10,a)") "@ECPMAIN: Unable to allocate",
     &         Mem_needed, " words of double precision memory"
        Call Errex
      Endif
C
      IECP_int      = Istart_Heap
      IECP_grad_int = IECP_int +      (1/2)*(NAOBASIS+1)*NAOBASIS

      If (Grads)
     &    Ilast  = IECP_grad_int + (3/2)*(NAOBASIS+1)*NAOBASIS*
     &                             Ncenters

      Call A3_ECP_MAIN(Iecp, Ncenters, Namat, Iqmstr, Jcostr,
     &                 Nucstr, Nrcstr, Jstrt, Jrs, Work(Ialpha), 
     &                 Nalpha, Work(Ipcoef), Npcoef, Cord, 
     &                 Charge, IGenby, Spherical, Decp_ints(IECP_int))

      If (Grads) Call A3_ECP_GRDINT_MAIN(Iecp, Ncenters, Namat,
     &                                   Iqmstr, Jcostr,
     &                                   Nucstr, Nrcstr, Jstrt, 
     &                                   Jrs, Work(Ialpha),
     &                                   Nalpha, Work(Ipcoef), 
     &                                   Npcoef, Cord, Charge,
     &                                   IGenby, Spherical,
     &                                   Decp_ints(IECP_grad_int),
     &                                   Naobasis)

      Return  
      End
