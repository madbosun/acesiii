



c#define _ADD_H_AT_THE_VERY_END
c#define _PRINTCRAP

c THIS ROUTINE LOADS THE AO INTEGRALS FROM THE CORRESPONDING
c INTEGRAL FILE AND CONSTRUCTS THE FOCK MATRIX (ARRAY F) USING
c THE DENSITY MATRIX SUPPLIED IN D (UHF VERSION).
c
c COLOUMB TERMS:
c
c    FA(I,J) = FA(I,J) + DT(K,L) (IJ|KL)
c    FA(K,L) = FA(K,L) + DT(I,J) (IJ|KL)
c
c    FB(I,J) = FB(I,J) + DT(K,L) (IJ|KL)
c    FB(K,L) = FB(K,L) + DT(I,J) (IJ|KL)
c
c EXCHANGE TERMS:
c
c    FA(I,K) = FA(I,K) - 1/2 DA(J,L) (IJ|KL)
c    FA(J,L) = FA(J,L) - 1/2 DA(I,K) (IJ|KL)
c    FA(I,L) = FA(I,L) - 1/2 DA(J,K) (IJ|KL)
c    FA(J,K) = FA(J,K) - 1/2 DA(I,L) (IJ|KL)
c
c    FB(I,K) = FB(I,K) - 1/2 DB(J,L) (IJ|KL)
c    FB(J,L) = FB(J,L) - 1/2 DB(I,K) (IJ|KL)
c    FB(I,L) = FB(I,L) - 1/2 DB(J,K) (IJ|KL)
c    FB(J,K) = FB(J,K) - 1/2 DB(I,L) (IJ|KL)
c
c  J. GAUSS AND J.F. STANTON, QTP 1993
c  A. YAU AND A. PERERA, 2000 - MOLCAS integral file support

      SUBROUTINE MKUHFF(FA,FB,DA,DB,DT,H,BUF,IBUF,NTOTAL,
     &                  NBAST,NBAS,IMAP,ILNBUF,LUINT,ADDH,
     &                  naobasfn,iuhf,scfks,scfksexact,
     &                  scfkslastiter,
     &                  V,z1,ksa,ksb,
     &                  screxc,scr,scr2,
     &                  valao,valgradao,totwt,
     &                  max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)

c   FA - SYMMETRY BLOCKED FOCK MATRIX (ALPHA)
c   FB - SYMMETRY BLOCKED FOCK MATRIX (BETA)
c   DA - SYMMETRY BLOCKED DENSITY MATRIX (ALPHA)
c   DB - SYMMETRY BLOCKED DENSITY MATRIX (BETA)
c   H  - SYMMETRY BLOCKED ONE-ELECTRON INTEGRALS
c
c   SEWARD:
c      BUF  - square-packed, unscaled density matrix
c             (vs. scaled, triangular-packed in D)
c      IBUF - work buffer for doubles, used like the array Work in
c             MOLCAS/ftwo.f
c             BEWARE: The compiler addresses this array as integers, so
c                     all the offset counters must be scaled by iintfp.
c                     In addition, you cannot (meaningfully) address any
c                     one element in the array, you must pass the
c                     address to a subroutine.
c      IMAP - NOT USED!!!
c      ILNBUF - maximum number of doubles that will fit in IBUF
c               (plug-in for GetMem(...,'MAX',...))
c   .NOT.SEWARD:
c      BUF  - BUFFER OF SIZE ILNBUF*IINTFP FOR INTEGRALS
c      IBUF - BUFFER OF SIZE ILNBUF FOR INDICES
c      IMAP - SYMMETRY VECTOR (FULL ARRAY -> SYMMETRY BLOCKED)
c      ILNBUF - BUFFER LENGTH
c
c   NTOTAL  - TOTAL LENGTH OF D AND F
c   NBAST   - TOTAL NUMBER OF BASIS FUNCTIONS
c   NBAS(I) - NUMBER OF BASIS FUNCTIONS IN IRREP I
c
c   LUINT - UNIT NUMBER FOR INTEGRAL FILE
c
c   ADDH - LOGICAL SWITCH FOR ADDING H TO F

c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
c
c ARGUMENT LIST
      INTEGER NTOTAL, NBAST, NBAS(8), ILNBUF, LUINT
      DOUBLE PRECISION FA(NTOTAL), FB(NTOTAL)
      DOUBLE PRECISION DA(NTOTAL), DB(NTOTAL), DT(NTOTAL)
      DOUBLE PRECISION H(NTOTAL), BUF(ILNBUF)
      INTEGER IBUF(ILNBUF), IMAP(NBAST,NBAST)
      LOGICAL ADDH
c
c---------------------------------------------------------------







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c This contains flags that are set in the INTGRT namelist.  See the
c file initintgrt.F for a description of each of them.
c
c The following are exceptions:
c    int_ks           : .true. if we are doing Kohn-Sham
c    int_ks_finaliter : .true. if this is the final iteration of KS
c    int_ks_exch      : which potential to use to calculate exchange
c    int_ks_corr      : which potential to use to calculate correlation
c    int_kspot        : which hybrid functional to use to calculate potential
c                       (if equal to fun_special, use int_ks_exch and
c                       int_ks_corr)
c    int_dft_fun      : which functional to use with any SCF density
c                       Added and modified by Stan Ivanov
c                       (if fun_special the functional is user-defined
c                        if fun_hyb_name then use hybrid functional) 
c    int_printlev     : 0 if we are doing a dft calculation, 1 if we
c                       are doing the final iteration of a KS calculation,
c                       2 if we are doing a KS iteration.
c These are set in the calling routines, NOT in the namelist.
c                     : Additions by S. Ivanov
c     num_acc_ks      : .true.  if numerical accelerator is used for KS 
c                        Default is .true.
c     ks_exact_ex     : .true. if exact LOCAL exchange is used for KS
c                        Deafult is .false.
c     int_tdks        : .true. if time-dependent KS calculation is
c                        requested
c                        Default is .false.
c     int_ks_scf      : .true. if the actual KS SCF energy is being
c                        calculated and printed out. Default is .false.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defenegrid,int_defenetype,
     &    int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp

      logical int_ks,num_acc_ks,ks_exact_ex,int_tdks,int_ks_scf,
     &        int_ks_finaliter

      double precision
     &    int_radlimit,coef_pot_nonlocal

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defenegrid,
     &    int_defenetype,int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_ks_finaliter,int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp
c

c  prakash added int_ksmem to the common block
      common /intgrtflagsd/ int_radlimit,coef_pot_nonlocal
      common /intgrtflagsl/ int_ks,num_acc_ks,ks_exact_ex,int_tdks,
     &                      int_ks_scf

      save /intgrtflags/
      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_prt_never,int_prt_dft,int_prt_ks,int_prt_always
      parameter (int_prt_never   =1)
      parameter (int_prt_dft     =2)
      parameter (int_prt_ks      =3)
      parameter (int_prt_always  =4)

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)


c
       integer iuhf,naobasfn,ncount,natoms,max_angpts,
     &         intnumradpts
       logical scfks,scfksexact,scfkslastiter,kshf
       integer z1(naobasfn,2)
       double precision V(naobasfn,naobasfn,iuhf+1),   
     &                  ksa(NTOTAL),ksb(ntotal),
     &                  screxc(naobasfn,naobasfn),
     &                  scr(naobasfn,naobasfn,iuhf+1),
     &                  scr2(nbast,nbast),coef_nonloc,
     &                  valao(naobasfn,intnumradpts,
     &                        max_angpts,ncount),
     &                  valgradao(naobasfn,intnumradpts,
     &                            max_angpts,ncount,3),
     &                  totwt(ncount,intnumradpts,max_angpts)
c------------------------------------------------------------------
c INTERNAL VARIABLES FOR/FROM VMOL
      DOUBLE PRECISION X
      INTEGER I, J, K, L, IOFF, IOFFT
      INTEGER IJ, IK, IL, JK, JL, KL, INDI, INDJ, INDK, INDL
      INTEGER IRREP, INDS, ILENGTH
      INTEGER NAOBUF, NUMINT, NUT
      CHARACTER*80 FNAME
c INTERNAL VARIABLES FOR/FROM MOLCAS
      double precision DDot, Temp
      integer nSym, BetaOff,
     &        iSym, kSym, nBs, iBs, jBs, kBs, mBs,
     &        iRc, iOpt, lbuf, nSub,
     &        ijPairs, klPairs, nInts, lw2,
     &        ipInts, ipTri, ipSqr, jpTri, jpSqr, kpTri, kpSqr

c PARAMETERS
      DOUBLE PRECISION TWOM, ONEM, HALFM, ZERO, FOURTH, HALF, ONE, TWO
      PARAMETER (TWOM=-2.0D0, ONEM=-1.0D0)
      PARAMETER (HALFM=-0.5D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (FOURTH=0.25D0, HALF=0.5D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)

c COMMON BLOCKS
c molcas.com : begin
      logical seward, petite_list
      character*8 fnord
      integer luord
      integer ipmat(8,2), lbbt, lbbs
      common /molcas_com/ seward, petite_list, fnord, luord,
     &                    ipmat, lbbt, lbbs
c molcas.com : end


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



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
cwc1


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





      common /GMSdirect/ dirscf,fdiff,schwrz,pople,i011,i012,i013
      logical dirscf,fdiff,schwrz,pople
      integer i011,i012,i013
cwc0
c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end
      INTEGER ITER
      COMMON /LOGBLOCK/ ITER

c FUNCTION DEFINITIONS
      INTEGER IAND,IOR,ISHFT
      INTEGER INT, IUPKI, IUPKJ, IUPKK, IUPKL, IPACK, INDX
      IUPKI(INT)=IAND(INT,IALONE)
      IUPKJ(INT)=IAND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=IAND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=IAND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=IOR(IOR(IOR(I,ISHFT(J,IBITWD)),
     &                             ISHFT(K,2*IBITWD)),
     &                             ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2

c---------------------------------------------------------------------
c prakash 23/04
           if ( scfks) then
               call dzero(ksa,ntotal)
               call dzero(ksb,ntotal)
            end if
c ----------------------------------------------------------------------

      if (petite_list) then
         nSym = NIRREP
c      o MOLCAS integrals are stored in upper-triangular, packed
c        form. Lucky for us, we too store the density and build the
c        Fock matrices in upper-packed form. A few notes:
c         1) At the end of VMOL integral processing, we scale the
c            whole Fock matrix up by 2 and then the on-diagonal
c            elements down by 1/2. Given the succinct way Seward
c            stores its integrals, we can accomplish all this work in
c            one (or two) fell swoop(s).
c         2) The exchange contributions use square density blocks
c            while the coulomb contributions use triangular-packed
c            blocks. Watch out for this if you intend on changing
c            stuff around.
c         3) I threw out the threshold processing found in the Molcas
c            SCF. If someone ever decides to put it back, then look
c            at $MOLCAS/src/scf/ftwo.f as a reference.
c      o In order to accomodate both Alpha and Beta square-packed
c        densities, we need to put them one-after-the-other in Buf().
c        The BetaOff variable USED WITH ANOTHER POINTER will accurately
c        point to elements in the square Beta density.
         BetaOff = ISQRLN(NIRREP+1)
         Do iSym = 1,nSym
            nBs = nBas(iSym)
            If (nBs.ne.0) Then
               ijPairs = ishft((nBs*(nBs+1)),-1)
c            o ALPHA
               Do i = 0,(ijPairs-1)
                  Buf(ipMat(iSym,lBBS)+i)=DA(ipMat(iSym,lBBT)+i)
               End Do
               Call SP2Sy(nBs,Buf(ipMat(iSym,lBBS)),nBs,'U',iOff)
               If (iOff.ne.0) Then
                  write(*,*) '@MKUHFF: There was a problem unpacking ',
     &                       'the alpha density.'
                  Call ErrEx
                  Stop 1
               End If
c            o BETA
               Do i = 0,(ijPairs-1)
                  Buf(BetaOff+ipMat(iSym,lBBS)+i)=DB(ipMat(iSym,lBBT)+i)
               End Do
               Call SP2Sy(nBs,Buf(BetaOff+ipMat(iSym,lBBS)),
     &                    nBs,'U',iOff)
               If (iOff.ne.0) Then
                  write(*,*) '@MKUHFF: There was a problem unpacking ',
     &                       'the beta density.'
                  Call ErrEx
                  Stop 1
               End If
            End If
         End Do
c      o collect DA and DB into DT
         Do i=1,ntotal
            DT(i)=DA(i)+DB(i)
         End Do
c      o initialize the Fock matrix
         If (addh) Then
            Do i=1,ntotal
               FA(i)=H(i)
            End Do
            Do i=1,ntotal
               FB(i)=H(i)
            End Do
         Else
            Do i=1,ntotal
               FA(i)=ZERO
            End Do
            Do i=1,ntotal
               FB(i)=ZERO
            End Do
         End If
      else
c      o initialize IMAP
         CALL IZERO(IMAP,NBAST*NBAST)
         IOFF=0
         IOFFT=0
         DO IRREP=1,NIRREP
            DO I=1,NBAS(IRREP)
               DO J=1,NBAS(IRREP)
                  INDS=INDX(MAX(I,J),MIN(I,J))+IOFFT
                  IMAP(I+IOFF,J+IOFF)=INDS
               END DO
            END DO
            IOFF=IOFF+NBAS(IRREP)
            IOFFT=IOFFT+NBAS(IRREP)*(NBAS(IRREP)+1)/2
         END DO
c      o initialize the Fock matrix
         DO I = 1,NTOTAL
            FA(I)=ZERO
         END DO
         DO I = 1,NTOTAL
            FB(I)=ZERO
         END DO
         CALL VADD(DT,DA,DB,NTOTAL,ONE)
cwc1
         if (dirscf) call SCFINT (.false.,FA,FB,DA,DB,icore(i013),
     >                 imap,nbast,ntotal,icore(i011),icore(i012))
cwc0
c     end if (seward)
      end if
c ----------------------------------------------------------------------

      if (petite_list) then

ccccccccccccccccccccccccccc
c START SEWARD PROCESSING c
ccccccccccccccccccccccccccc

c   o This stuff is ONLY for Coulomb contributions. It is unclear whether
c     the calling program uses the total density on output, so we will
c     unscale the array at the end to be safe.
      Do iSym = 1,nSym
         Call dDiagScal_SP('U','O',nBas(iSym),TWO,DT(ipMat(iSym,lBBT)))
      End Do

c   o (ii|ii) integrals
      Do iSym=1,nSym
         nBs=nBas(iSym)
         ijPairs=nBs*(nBs+1)/2
         klPairs=ijPairs
         nInts=ijPairs*klPairs
         If ( nInts.ne.0 ) Then
            lw2=(nBs*nBs+1)*iintfp
            lBuf=MIN(iLnBuf,nInts+1)
            lBuf=MAX(lBuf,klPairs+1)
            iOpt=1
            Call RdOrd(iRc,iOpt,iSym,iSym,iSym,iSym,
     &                 iBuf(lw2),lBuf,nSub)
            iOpt=2
            ipInts=0
            ipTri=ipMat(iSym,lBBT)
            ipSqr=ipMat(iSym,lBBS)
            Do iBs=1,nBs
               jpTri=ipMat(iSym,lBBT)
               jpSqr=ipMat(iSym,lBBS)
               Do jBs=1,iBs
                  If ( nSub.eq.0 ) Then
                     Call RdOrd(iRc,iOpt,iSym,iSym,iSym,iSym,
     &                          iBuf(lw2),lBuf,nSub)
                     ipInts=0
                  End If
                  iOff=ipTri+jBs-1
                  Temp=DDot(klPairs,iBuf(lw2+ipInts),1,
     &                              DT(ipMat(iSym,lBBT)),1)
                  Call Square(iBuf(lw2+ipInts),iBuf(1),1,nBs,nBs)
c               o ALPHA
                  FA(iOff)=FA(iOff)+Temp
                  Call DGeMV('N',jBs,nBs,
     &                       ONEM,iBuf(1),nBs,
     &                            Buf(ipSqr),1,
     &                       ONE, FA(jpTri),1)
                  If ( iBs.ne.jBs ) Call
     &               DGeMV('N',iBs,nBs,
     &                     ONEM,iBuf(1),nBs,
     &                          Buf(jpSqr),1,
     &                     ONE, FA(ipTri),1)
c               o BETA
                  FB(iOff)=FB(iOff)+Temp
                  Call DGeMV('N',jBs,nBs,
     &                       ONEM,iBuf(1),nBs,
     &                            Buf(BetaOff+ipSqr),1,
     &                       ONE, FB(jpTri),1)
                  If ( iBs.ne.jBs ) Call
     &               DGeMV('N',iBs,nBs,
     &                     ONEM,iBuf(1),nBs,
     &                          Buf(BetaOff+jpSqr),1,
     &                     ONE, FB(ipTri),1)
                  nSub=nSub-1
                  ipInts=ipInts+klPairs*iintfp
                  jpTri=jpTri+jBs
                  jpSqr=jpSqr+nBs
c              End Do jBs=1,iBs
               End Do
               ipTri=ipTri+iBs
               ipSqr=ipSqr+nBs
c           End Do iBs=1,nBs
            End Do
c        End If ( nInts.ne.0 )
         End If
c     End Do iSym=1,nSym
      End Do

c   o (ii|jj) integrals, i>j
      Do iSym=2,nSym
         nBs=nBas(iSym)
         ijPairs=ishft((nBs*(nBs+1)),-1)
         Do kSym=1,iSym-1
            mBs=nBas(kSym)
            klPairs=ishft((mBs*(mBs+1)),-1)
            nInts=ijPairs*klPairs
            If ( nInts.ne.0 ) Then
               lBuf=MIN(iLnBuf,nInts+1)
               lBuf=MAX(lBuf,klPairs+1)
               iOpt=1
               Call RdOrd(iRc,iOpt,iSym,iSym,kSym,kSym,
     &                    iBuf(1),lBuf,nSub)
               iOpt=2
               ipInts=0
               ipTri=ipMat(iSym,lBBT)
               Do iBs=1,nBs
                  Do jBs=1,iBs
                     iOff=ipTri+jBs-1
                     If ( nSub.eq.0 ) Then
                        Call RdOrd(iRc,iOpt,iSym,iSym,kSym,kSym,
     &                             iBuf(1),lBuf,nSub)
                        ipInts=0
                     End If
                     Temp=DDot(klPairs,iBuf(1+ipInts),1,
     &                                 DT(ipMat(kSym,lBBT)),1)
                     FA(iOff)=FA(iOff)+Temp
                     FB(iOff)=FB(iOff)+Temp
                     Temp=DT(iOff)
                     Call dAXPY(klPairs,Temp,
     &                          iBuf(1+ipInts),1,
     &                          FA(ipMat(kSym,lBBT)),1)
                     Call dAXPY(klPairs,Temp,
     &                          iBuf(1+ipInts),1,
     &                          FB(ipMat(kSym,lBBT)),1)
                     nSub=nSub-1
                     ipInts=ipInts+klPairs*iintfp
                  End Do
                  ipTri=ipTri+iBs
               End Do
c           End If ( nInts.ne.0 )
            End If
c        End Do kSym=1,iSym-1
         End Do
c     End Do iSym=2,nSym
      End Do

c   o Now that the Coulomb contributions are done, we "fix" the density.
c     I am not sure if we even need to do this since DT is never used
c     again in vscf, but there is no statement (such as INTENT=OUT)
c     that would suggest DT's permanence.
      Do iSym = 1,nSym
         Call dDiagScal_SP('U','O',nBas(iSym),HALF,DT(ipMat(iSym,lBBT)))
      End Do

c   o (ij|ij) integrals, i>j
      Do iSym=2,nSym
         nBs=nBas(iSym)
         Do kSym=1,iSym-1
            mBs=nBas(kSym)
            ijPairs=nBs*mBs
            klPairs=ijPairs
            nInts=ijPairs*klPairs
            If ( nInts.ne.0 ) Then
               lBuf=MIN(iLnBuf,nInts+1)
               lBuf=MAX(lBuf,klPairs+1)
               iOpt=1
               Call RdOrd(iRc,iOpt,iSym,kSym,iSym,kSym,
     &                    iBuf(1),lBuf,nSub)
               iOpt=2
               ipInts=0
               ipSqr=ipMat(iSym,lBBS)
               ipTri=ipMat(iSym,lBBT)
               Do iBs=1,nBs
                  kpSqr=ipMat(kSym,lBBS)
                  kpTri=ipMat(kSym,lBBT)
                  Do kBs=1,mBs
                     If ( nSub.eq.0 ) Then
                        Call RdOrd(iRc,iOpt,iSym,kSym,iSym,kSym,
     &                             iBuf(1),lBuf,nSub)
                        ipInts=0
                     End If
c                  o ALPHA
                     Call DGeMV('N',kBs,nBs,
     &                          ONEM,iBuf(1+ipInts),mBs,
     &                               Buf(ipSqr),1,
     &                          ONE, FA(kpTri),1)
                     Call DGeMV('T',mBs,iBs,
     &                          ONEM,iBuf(1+ipInts),mBs,
     &                               Buf(kpSqr),1,
     &                          ONE, FA(ipTri),1)
c                  o BETA
                     Call DGeMV('N',kBs,nBs,
     &                          ONEM,iBuf(1+ipInts),mBs,
     &                               Buf(BetaOff+ipSqr),1,
     &                          ONE, FB(kpTri),1)
                     Call DGeMV('T',mBs,iBs,
     &                          ONEM,iBuf(1+ipInts),mBs,
     &                               Buf(BetaOff+kpSqr),1,
     &                          ONE, FB(ipTri),1)
                     nSub=nSub-1
                     ipInts=ipInts+klPairs*iintfp
                     kpSqr=kpSqr+mBs
                     kpTri=kpTri+kBs
                  End Do
                  ipSqr=ipSqr+nBs
                  ipTri=ipTri+iBs
               End Do
c           End If ( nInts.ne.0 )
            End If
c        End Do kSym=1,iSym-1
         End Do
c     End Do iSym=2,nSym
      End Do

c   o This is only for collecting terms from a parallel process.
c      Call GAdsum(F,NTOTAL)

ccccccccccccccccccccccccc
c END SEWARD PROCESSING c
ccccccccccccccccccccccccc

      else
cwc1
         if (dirscf) goto 999
cwc0

ccccccccccccccccccccccccc
c START VMOL PROCESSING c
ccccccccccccccccccccccccc

c IIII CONTRIBUTES TO BOTH COLOUMB AND EXCHANGE INTEGRALS
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      NAOBUF=0
      NUMINT=0
      CALL LOCATE(LUINT,'TWOELSUP')
 1    CONTINUE
         READ(LUINT) BUF, IBUF, NUT
         NAOBUF=NAOBUF+1
         DO INT=1,NUT
            X=BUF(INT)
            INDI=IUPKI(IBUF(INT))
            INDJ=IUPKJ(IBUF(INT))
            INDK=IUPKK(IBUF(INT))
            INDL=IUPKL(IBUF(INT))
            IJ=IMAP(INDJ,INDI)
            KL=IMAP(INDL,INDK)
            IK=IMAP(INDK,INDI)
            JL=IMAP(INDL,INDJ)
            IL=IMAP(INDL,INDI)
            JK=IMAP(INDJ,INDK)
            IF (INDI.EQ.INDJ) X=X*HALF
            IF (INDK.EQ.INDL) X=X*HALF
            IF (IJ.EQ.KL)     X=X*HALF
            FA(IJ)=FA(IJ)+DT(KL)*X
            FA(KL)=FA(KL)+DT(IJ)*X
            FB(IJ)=FB(IJ)+DT(KL)*X
            FB(KL)=FB(KL)+DT(IJ)*X

             if (scfks.AND. .not.scfkslastiter) then
                if(kshf) then
                   coef_pot_nonlocal=0.d0
                end if
                 if(coef_pot_nonlocal .ne.0.d0)then
                    coef_nonloc=coef_pot_nonlocal
                 else
                    goto 50
                 end if
             else
                 coef_nonloc=1.D0
            end if
 
            FA(IK)=FA(IK)+coef_nonloc*HALFM*DA(JL)*X
            FA(JL)=FA(JL)+coef_nonloc*HALFM*DA(IK)*X
            FA(IL)=FA(IL)+coef_nonloc*HALfM*DA(JK)*X
            FA(JK)=FA(JK)+coef_nonloc*HALFM*DA(IL)*X
            FB(IK)=FB(IK)+coef_nonloc*HALFM*DB(JL)*X
            FB(JL)=FB(JL)+coef_nonloc*HALFM*DB(IK)*X
            FB(IL)=FB(IL)+coef_nonloc*HALfM*DB(JK)*X
            FB(JK)=FB(JK)+coef_nonloc*HALFM*DB(IL)*X

  50     continue

         END DO
         NUMINT=NUMINT+NUT
      IF (NUT.NE.-1) GOTO 1
      CLOSE(UNIT=LUINT,STATUS='KEEP')

c FOR NIRREP > 1, READ ALSO IIJJ AND IJIJ INTEGRAL FILES
      IF (NIRREP.GT.1) THEN

c IIJJ CONTRIBUTES ONLY TO COLOUMB INTEGRALS
         CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
         OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &        ACCESS='SEQUENTIAL')
         NAOBUF=0
         NUMINT=0
         CALL LOCATE(LUINT,'TWOELSUP')
 101     CONTINUE
            READ(LUINT) BUF, IBUF, NUT
            NAOBUF=NAOBUF+1
            DO INT=1,NUT
               X=BUF(INT)
               INDI=IUPKI(IBUF(INT))
               INDJ=IUPKJ(IBUF(INT))
               INDK=IUPKK(IBUF(INT))
               INDL=IUPKL(IBUF(INT))
               IJ=IMAP(INDJ,INDI)
               KL=IMAP(INDL,INDK)
               IF (INDI.EQ.INDJ) X=X*HALF
               IF (INDK.EQ.INDL) X=X*HALF
               FA(IJ)=FA(IJ)+DT(KL)*X
               FA(KL)=FA(KL)+DT(IJ)*X
               FB(IJ)=FB(IJ)+DT(KL)*X
               FB(KL)=FB(KL)+DT(IJ)*X
            END DO
            NUMINT=NUMINT+NUT
         IF (NUT.NE.-1) GOTO 101
         CLOSE(UNIT=LUINT,STATUS='KEEP')

c IJIJ CONTRIBUTES ONLY TO EXCHANGE
         CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
         OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &        ACCESS='SEQUENTIAL')
         NAOBUF=0
         NUMINT=0
         CALL LOCATE(LUINT,'TWOELSUP')
 201     CONTINUE
            READ(LUINT) BUF, IBUF, NUT
            NAOBUF=NAOBUF+1
            DO INT=1,NUT
               X=BUF(INT)
               INDI=IUPKI(IBUF(INT))
               INDJ=IUPKJ(IBUF(INT))
               INDK=IUPKK(IBUF(INT))
               INDL=IUPKL(IBUF(INT))
               IJ=INDX(MAX(INDJ,INDI),MIN(INDJ,INDI))
               KL=INDX(MAX(INDL,INDK),MIN(INDL,INDK))
               IK=IMAP(INDK,INDI)
               JL=IMAP(INDL,INDJ)
               IF (IJ.EQ.KL) X=X*HALF

C____________________________________________________________________
                if (scfks  .and. .not. scfkslastiter) then
                   if(kshf) then
                     coef_pot_nonlocal=0.d0
                   end if           
                   if (coef_pot_nonlocal .ne.0d0) then
                       coef_nonloc=coef_pot_nonlocal
                   else
                       goto 60
                   end if
                else
                   coef_nonloc=1.d0
                end if

               FA(IK)=FA(IK)+coef_nonloc*HALFM*DA(JL)*X
               FA(JL)=FA(JL)+coef_nonloc*HALFM*DA(IK)*X
               FB(IK)=FB(IK)+coef_nonloc*HALFM*DB(JL)*X
               FB(JL)=FB(JL)+coef_nonloc*HALFM*DB(IK)*X

   60        continue
 
            END DO
            NUMINT=NUMINT+NUT
         IF (NUT.NE.-1) GOTO 201
         CLOSE(UNIT=LUINT,STATUS='KEEP')

c     END IF (NIRREP.GT.1)
      END IF
cwc1
 999  continue

cwc0
c SCALE DIAGONAL PARTS OF F BY A FACTOR OF TWO
      CALL SCALEF(FA,NTOTAL,NBAS)
      CALL SCALEF(FB,NTOTAL,NBAS)
      CALL XSCAL(NTOTAL,TWO,FA,1)
      CALL XSCAL(NTOTAL,TWO,FB,1)


C-----------------------------------------------------------------
C ks terms
         if(scfks.AND. .not. scfkslastiter) then
         if(iter.gt.1)then
           call dzero(V,naobasfn*naobasfn*(1+iuhf))
           call integxc(V,z1,screxc,scr,scfkslastiter,
     &                valao,valgradao,totwt,
     &       intnumradpts,max_angpts,ncount,kshf)
           call mat_trans(0,2,v(1,1,1),ksa,0)
           call mat_trans(0,2,v(1,1,2),ksb,0)
           call daxpy(NTOTAL,one,ksa,1,FA,1)
           call daxpy(NTOTAL,one,ksb,1,FB,1)
         end if
         end if
C end prakash
C________________________________________________________________

c ADD THE 1-ELECTRON CONTRIBUTION
      IF (ADDH) THEN
         CALL XAXPY(NTOTAL,ONE,H,1,FA,1)
         CALL XAXPY(NTOTAL,ONE,H,1,FB,1)
      END IF

ccccccccccccccccccccccc
c END VMOL PROCESSING c
ccccccccccccccccccccccc

c     end if (seward)
      end if


c      if (ITER.eq.2) stop 1

c
      RETURN
      END

