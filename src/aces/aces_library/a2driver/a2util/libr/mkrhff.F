c THIS ROUTINE LOADS THE AO INTEGRALS FROM THE CORRESPONDING
c INTEGRAL FILE AND CONSTRUCTS THE FOCK MATRIX (ARRAY F) USING
c THE DENSITY MATRIX SUPPLIED IN D.
c
c COLOUMB TERMS:
c
c    F(I,J) = F(I,J) + D(K,L) (IJ|KL)
c    F(K,L) = F(K,L) + D(I,J) (IJ|KL)
c
c EXCHANGE TERMS:
c
c    F(I,K) = F(I,K) - 1/4 D(J,L) (IJ|KL)
c    F(J,L) = F(J,L) - 1/4 D(I,K) (IJ|KL)
c    F(I,L) = F(I,L) - 1/4 D(J,K) (IJ|KL)
c    F(J,K) = F(J,K) - 1/4 D(I,L) (IJ|KL)
c
c  J. GAUSS AND J.F. STANTON, QTP 1993
c  A. YAU AND A. PERERA, 2000 - MOLCAS integral file support

      SUBROUTINE MKRHFF(F,D,H,BUF,IBUF,NTOTAL,NBAST,NBAS,
     &                  IMAP,ILNBUF,LUINT,ADDH,
     &                  naobasfn,iuhf,scfks,scfksexact,
     &                  scfkslastiter,iter,
     &                  V,z1,ksa,
     &                  screxc,scr,scr2,
     &                  screval,valao,valgradao,totwt,
     &                  max_angpts,natoms,intnumradpts,ncount,kshf)

c   F - SYMMETRY BLOCKED FOCK MATRIX
c   D - SYMMETRY BLOCKED DENSITY MATRIX
c   H - SYMMETRY BLOCKED ONE-ELECTRON INTEGRALS
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
C
C-----------------------------------------------------------------     
c ARGUMENT LIST
      INTEGER NTOTAL, NBAST, NBAS(8), ILNBUF, LUINT
      DOUBLE PRECISION F(NTOTAL), D(NTOTAL), H(NTOTAL), BUF(ILNBUF)
      INTEGER IBUF(ILNBUF), IMAP(NBAST,NBAST)
      LOGICAL ADDH

C--------------------------------------------------------------
C KS related decalartions, prakash 23/04
#include<aces.h>
#include<intgrtflags.com>
C
       integer iuhf,iter,naobasfn,ncount,natoms,max_angpts,
     &         intnumradpts
       logical scfks,scfksexact,scfkslastiter,kshf
C
       integer z1(naobasfn,2)
C
       double precision V(naobasfn,naobasfn,iuhf+1),
     &                  ksa(NTOTAL),
     &                  screxc(naobasfn,naobasfn),
     &                  scr(naobasfn,naobasfn,iuhf+1),
     &                  scr2(nbast,nbast),
     &                  screval(nbast),coef_nonloc, 
     &                  valao(naobasfn,max_angpts,
     &                        intnumradpts,ncount),
     &                  valgradao(naobasfn,max_angpts,
     &                            intnumradpts,ncount,3),
     &                  totwt(ncount,intnumradpts,max_angpts)
C-----------------------------------------------------------------
c INTERNAL VARIABLES FOR/FROM VMOL
      DOUBLE PRECISION X
      INTEGER I, J, K, L, IOFF, IOFFT
      INTEGER IJ, IK, IL, JK, JL, KL, INDI, INDJ, INDK, INDL
      INTEGER IRREP, INDS, ILENGTH
      INTEGER NAOBUF, NUMINT, NUT
      CHARACTER*80 FNAME
c INTERNAL VARIABLES FOR/FROM MOLCAS
      double precision DDot, Temp
      integer nSym,
     &        iSym, kSym, nBs, iBs, jBs, kBs, mBs,
     &        iRc, iOpt, lbuf, nSub,
     &        ijPairs, klPairs, nInts, lw2,
     &        ipInts, ipTri, ipSqr, jpTri, jpSqr, kpTri, kpSqr

c COMMON BLOCKS
#include "molcas.com"
#include "machsp.com"
#include "syminf.com"
#include "flags.com"
cwc1
#include "icore.com"
      common /GMSdirect/ dirscf,fdiff,schwrz,pople,i011,i012,i013
      logical dirscf,fdiff,schwrz,pople
      integer i011,i012,i013
cwc0

c PARAMETERS
      DOUBLE PRECISION HALFM, FOURTHM, ZERO, HALF, ONE, TWO
      PARAMETER (HALFM=-0.5D0, FOURTHM=-0.25D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)

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
C
C---------------------------------------------------------------------
c prakash 23/04
           if ( scfks) then
               call dzero(ksa,ntotal)
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
         Do iSym = 1,nSym
            nBs = nBas(iSym)
            If (nBs.ne.0) Then
               ijPairs = ishft((nBs*(nBs+1)),-1)
               Do i = 0,(ijPairs-1)
                  Buf(ipMat(iSym,lBBS)+i)=D(ipMat(iSym,lBBT)+i)
               End Do
               Call SP2Sy(nBs,Buf(ipMat(iSym,lBBS)),nBs,'U',iOff)
               If (iOff.ne.0) Then
                  write(*,*) '@MKRHFF: There was a problem unpacking ',
     &                       'the density.'
                  Call ErrEx
                  Stop 1
               End If
            End If
         End Do
c      o initialize the Fock matrix
         If (addh) Then
            Do i=1,ntotal
               F(i)=H(i)
            End Do
         Else
            Do i=1,ntotal
               F(i)=ZERO
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
         DO I=1,NTOTAL
            F(I)=ZERO
         END DO
cwc1
         if (dirscf) call SCFINT (.true.,F,F,D,D,icore(i013),
     >            imap,nbast,ntotal,icore(i011),icore(i012))
cwc0
      end if

c ----------------------------------------------------------------------

      if (petite_list) then
c      print *, '@MKRHFF: Using the ORDINT petite list.'

ccccccccccccccccccccccccccc
c START SEWARD PROCESSING c
ccccccccccccccccccccccccccc

c   o This stuff is ONLY for Coulomb contributions. We will undo it
c     right after the (II|JJ) block.
      Do iSym = 1,nSym
         Call dDiagScal_SP('U','O',nBas(iSym),TWO,D(ipMat(iSym,lBBT)))
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
                  F(iOff)=F(iOff)
     &                   +DDot(klPairs,iBuf(lw2+ipInts),1,
     &                                 D(ipMat(iSym,lBBT)),1)
                  Call Square(iBuf(lw2+ipInts),iBuf(1),1,nBs,nBs)
                  Call DGeMV('N',jBs,nBs,
     &                       HALFM,iBuf(1),nBs,
     &                             Buf(ipSqr),1,
     &                       ONE,  F(jpTri),1)
                  If ( iBs.ne.jBs )
     &               Call DGeMV('N',iBs,nBs,
     &                          HALFM,iBuf(1),nBs,
     &                                Buf(jpSqr),1,
     &                          ONE,  F(ipTri),1)
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
                     F(iOff)=F(iOff)
     &                      +DDot(klPairs,iBuf(1+ipInts),1,
     &                                    D(ipMat(kSym,lBBT)),1)
                     Temp=D(iOff)
                     Call dAXPY(klPairs,Temp,
     &                          iBuf(1+ipInts),1,
     &                          F(ipMat(kSym,lBBT)),1)
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
      Do iSym = 1,nSym
         Call dDiagScal_SP('U','O',nBas(iSym),HALF,D(ipMat(iSym,lBBT)))
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
                     Call DGeMV('N',kBs,nBs,
     &                          HALFM,iBuf(1+ipInts),mBs,
     &                                Buf(ipSqr),1,
     &                          ONE,  F(kpTri),1)
                     Call DGeMV('T',mBs,iBs,
     &                          HALFM,iBuf(1+ipInts),mBs,
     &                                Buf(kpSqr),1,
     &                          ONE,  F(ipTri),1)
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
            F(IJ)=F(IJ)+D(KL)*X
            F(KL)=F(KL)+D(IJ)*X

C prakash 23/04
C For B3LYP type DFT we need a paramter for non local exhange 
C coef_pot_nonlocal<=coef_pot_exch(fun_exch_hf)
C fun_exch_hf<=b3lypa
C b3lypa< dftfunc.com as a parameter
C
               if (scfks.AND..not.scfkslastiter) then
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
                                
            F(IK)=F(IK)+coef_nonloc*FOURTHM*D(JL)*X
            F(JL)=F(JL)+coef_nonloc*FOURTHM*D(IK)*X
            F(IL)=F(IL)+coef_nonloc*FOURTHM*D(JK)*X
            F(JK)=F(JK)+coef_nonloc*FOURTHM*D(IL)*X

  50     continue
C ____________________________________________________________________


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
               F(IJ)=F(IJ)+D(KL)*X
               F(KL)=F(KL)+D(IJ)*X
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
               IJ=INDX(INDJ,INDI)
               KL=INDX(INDL,INDK)
               IK=IMAP(INDK,INDI)
               JL=IMAP(INDL,INDJ)
               IF (IJ.EQ.KL) X=X*HALF

C____________________________________________________________________

                if (scfks.AND. .not. scfkslastiter) then
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
               F(IK)=F(IK)+coef_nonloc*FOURTHM*D(JL)*X
               F(JL)=F(JL)+coef_nonloc*FOURTHM*D(IK)*X

   60       continue
C___________________________________________________________________
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
      CALL SCALEF(F,NTOTAL,NBAS)
      CALL XSCAL(NTOTAL,TWO,F,1)
C-----------------------------------------------------------------
C ks terms
      if(scfks.AND. .not. scfkslastiter) then
         if(iter.gt.1)then
           call dzero(V,naobasfn*naobasfn*(1+iuhf))
           call integxc(V,z1,screxc,scr,scfkslastiter,
     &                       valao,valgradao,totwt,
     &        intnumradpts,max_ANGPTS,ncount,kshf)
           call mat_trans(MAT_FULL,MAT_TRI,v(1,1,1),ksa,0)
           call B_AXPY(NTOTAL,one,ksa,1,F,1)
         end if
       end if
C end prakash
C________________________________________________________________

C 
c ADD THE 1-ELECTRON CONTRIBUTION
      IF (ADDH) CALL XAXPY(NTOTAL,ONE,H,1,F,1)
#ifdef _PRINTCRAP
      write(*,*) "The fock matrix:"
      call realprt(F,ntotal)
#endif /* _PRINTCRAP */

ccccccccccccccccccccccc
c END VMOL PROCESSING c
ccccccccccccccccccccccc

c     end if (seward)
      end if

      RETURN
      END

