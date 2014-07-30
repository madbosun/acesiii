#include "flags.h"
      SubrOutine projec_IFC(HESS, PHESS, Tmp, PMAT, GMATRX_N, 
     &                      GMATRX_M, GRD, GRDTMP, NXM6)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Logical Constrain_opt

C Built the projected Hessian and Gradients for redundent internal
C optimizations (see JCP, 117, 9160, 2002). 
C
#include <machsp.com>
#include <mxatms.par>
#include <jodaflags.com>
C
      Dimension HESS(NXM6, NXM6), PMAT(NXM6, NXM6), PHESS(NXM6, NXM6),
     &          GMATRX_M(NXM6, NXM6), GMATRX_N(NXM6, NXM6),
     &          GRD(NXM6), Scr(Maxredunco), Tmp(NXM6, NXm6),
     &          GRDTM(NXM6)
 
      LENGMAT=NXM6*NXM6
      CALL DGETREC(20,'JOBARC','GI-MATRX',LENGMAT,GMATRX_M)
      CALL DGETREC(20,'JOBARC','G-MATRX ',LENGMAT,GMATRX_N)

      Call XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, GMATRX_M, 
     &            NXM6, GMATRX_N, NXM6, 0.0D0, PMAT, NXM6)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The Hessian projector:P"
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
C Constrained redundent Internal optimizations; Built the C matrix
C as described in JCC, 17, 49, 1996. This applies to constrained 
C redundent internal optimiztions. Use Phes to store the C matrix,
C 
      Call IGetrec(0, "JOBARC", "CONSTRNS", Ilength, Ijunk)
      Constrain_opt = .FALSE.
C
#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,I4,1x,I4)") "Constrained opt?", Ilength,
     &                            Iflags2(h_IFLAGS2_opt_control)
#endif
C
      If (Ilength .Gt. 0) Then
          Constrain_opt = .TRUE.
          If (Iflags2(h_IFLAGS2_opt_control) .NE. 3) Then
             Write(6,"(3x,2a)") "Inconsitency: There are constrains in", 
     &            " input, but opt_control=constrained is not set."
              Call Errex
           Endif
      Endif
C 
      If (Constrain_opt) Then

         Call Build_cmatrix(NXM6, Phess)
C
C Built the CPC; Not that Gmatrx_N{I} are used as scratch matrices.
C     
         CALL XGEMM("T", "N", NXM6, NXM6, NXM6, 1.0D0, Phess,
     &               NXM6, Pmat, NXM6, 0.0D0, GMATRX_N, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, GMATRX_N,
     &               NXM6, Phess, NXM6, 0.0D0, GMATRX_M, NXM6)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The CPC"
      CALL OUTPUT(GMATRX_M, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
C Invert the CPC 
C
CSSS         Call Dgefa(GMATRX_M, NXM6, NXM6, Ipvt, Info)
C
C#ifdef _NOSKIP
CSSS      CALL MINV(GMATRX_M, NXM6, NXM6, Tmp, Det, 1.0D-8, 0, 1)
      
      CALL EIG(GMATRX_M, Tmp, NXM6, NXM6, 1)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The eigen vectors of the CPC"
      CALL OUTPUT(Tmp, 1, nxm6, 1, nxM6, nxm6, nxm6, 1)
      Write(6,"(a)") "The eigenvalues of the CPC "
      Write(6, "(4F10.5)") (Gmatrx_M(I,I), I=1, Nxm6)
#endif

      DO I = 1, NXM6
         IF (GMATRX_M(I, I) .LE. 1.0D-10) THEN
             GMATRX_M(I, I) = 0.0D0
         ELSE
             GMATRX_M(I, I) = 1.0D0/GMATRX_M(I, I)
         ENDIF
      ENDDO
C
C Built the generalized inverse of CPC matrix,
C
      CALL XGEMM('N', 'N', NXM6, NXM6, NXM6, 1.0D0, Tmp,
     &            NXM6, GMATRX_M, NXM6, 0.0D0, GMATRX_N, NXM6)
      CALL XGEMM('N', 'T', NXM6, NXM6, NXM6, 1.0D0, GMATRX_N,
     &            NXM6, Tmp, NXM6, 0.0D0, GMATRX_M, NXM6)
C#endif 
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The inverse of CPC"
      CALL OUTPUT(GMATRX_M, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
C Create CP and PC
C
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Phess,
     &               NXM6, Pmat, NXM6, 0.0D0, GMATRX_N, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Pmat,
     &               NXM6, Phess, NXM6, 0.0D0, Tmp, NXM6)
C
C Built PC(CPC)^(-1)CP
C 
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Tmp,
     &               NXM6, GMATRX_M, NXM6, 0.0D0, Phess, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Phess,
     &               NXM6, GMATRX_N, NXM6, 0.0D0, GMATRX_M, NXM6)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The PC(CPC)^(-1)CP"
      CALL OUTPUT(GMATRX_M, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
C Built P = P - PC(CPC)^(-1)CP

       Call XDaxpy(NXM6*NXM6, -1.0D0, GMATRX_M, 1, Pmat, 1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The constrained Hessian projector:P"
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
C Built PHP (projected Hessian)
C
     
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)

         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)

C Also project the gradients ie. form PG. Then copy the projected gradient
C to the FI array so that the rest of the code can proceed.
C
      Write(6,"(a)") "The gradients, g"
      Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NXM6)
         CALL XGEMM('N', 'N', NXM6, 1, NXM6, 1.0D0, PMAT,
     &               NXM6, GRD, NXM6, 0.0D0, GRDTMP, NXM6)
         CALL XDCOPY(NXM6, GRDTMP, 1, GRD, 1)
C
#ifdef _NOSKIP
C
C Built Alpha(1-P) where Alpha = 1000.0
      CALL ZERO(TMP,NXM6*NXM6)
      DO I=1, NXM6
         TMP(I,I) = 1.0d0
      END DO
C
C Let's do the PHP + 1000.0(1-P) (note that we have already
C built the PHP and stored in HESS). Make -(1-P)
C
      CALL XDAXPY(NXM6*NXM6,-1.d0,TMP,1,PMAT,1)
      CALL XDSCAL(NXM6*NXM6,1000.0D0,PMAT,1)
      CALL XDAXPY(NXM6*NXM6,-1.d0,PMAT,1,HESS,1)
#endif
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The projected Hessian, (PHP)"
      CALL OUTPUT(HESS, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
      Write(6,"(a)") "The projected gradients, Hg"
      Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NXM6)
#endif
      
      Else 
C
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)
C
C Also project the gradients. Then copy the projected gradient
C to the FI array so that the rest of the code can proceed.
C
         CALL XGEMM('N', 'N', NXM6, 1, NXM6, 1.0D0, PMAT,
     &               NXM6, GRD, NXM6, 0.0D0, GRDTMP, NXM6)
         CALL XDCOPY(NXM6, GRDTMP, 1, GRD, 1)

#ifdef _NOSKIP
      CALL ZERO(TMP,NXM6*NXM6)
      DO I=1, NXM6
         TMP(I,I) = 1.0d0
      END DO
C
C Let's do the PHP + 1000.0(1-P) (note that we have already
C built the PHP and stored in HESS).
C
      CALL XDAXPY(NXM6,-1.d0,TMP,1,PMAT,1)
      CALL XDSCAL(NXM6*NXM6,1000.0D0,PMAT,1)
      CALL XDAXPY(NXM6*NXM6,-1.d0,PMAT,1,HESS,1)
C
#endif
C
#ifdef _DEBUG_LVL0
          Write(6,"(a)")"The projected hessian (PHP)"
          CALL OUTPUT(HESS, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
          Write(6,*) 
          Write(6,"(a)") "The projected gradients Pg"
          Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NXM6)
#endif
      Endif 
           
C------The lines below are for Debugging only--------------

#ifdef _DEBUG_HESSIAN_PROJC
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The (PHP)(PHP)t matrix"
      CALL OUTPUT(HESS, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
      CALL EIG(HESS, PHESS, NXM6, NXM6, 1)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,"(a)") "The eigen vectors of the PHP"
      CALL OUTPUT(phess, 1, nxm6, 1, nxM6, nxm6, nxm6, 1)
      Write(6,"(a)") "The eigenvalues of the PHP "
      Write(6, "(4F10.5)") (Hess(I,I), I=1, Nxm6)
#endif
      DO I = 1, NXM6
         IF (HESS(I, I) .LE. 1.0D-10) THEN
             HESS(I, I) = 0.0D0
         ELSE
             HESS(I, I) = 1.0D0/HESS(I, I)
         ENDIF
      ENDDO
C
C Built the generalized inverse of H-matrix,
C
      CALL XGEMM('N', 'N', NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, HESS, NXM6, 0.0D0, GMATRX_N, NXM6)
      CALL XGEMM('N', 'T', NXM6, NXM6, NXM6, 1.0D0, GMATRX_N,
     &            NXM6, PHESS, NXM6, 0.0D0, HESS, NXM6)
C
C Built P(PHP)-P 
C
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)
C
C---These steps require defining STEP and TMP arrays.
C
      CALL XGEMM("N", "N", NXM6, 1, NXM6, -1.0D0, HESS,
     &            NXM6, GRD, NXM6, 0.0D0, STEP, NXM6)
      CALL IGETREC(20,'JOBARC','NUMREDCO', 1, NULLEVAL)
      CALL DGETREC(20,'JOBARC','REDEVECS',NULLEVAL*NXM6,
     &            GMATRX_N)
      CALL XGEMM('T', 'N', NULLEVAL, 1, NXM6, 1.0D0,
     &               GMATRX_N, NXM6, STEP, 15,
     &               0.0D0, TMP1, 3)
     
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The step and overlap with the red. space"
      Write(6, "(8F10.5)") (Step(I), I=1, Nxm6)
      Write(6, "(8F10.5)") (TMP1(I), I=1, NULLEVAL)
#endif
      STOP
#endif
C------------------------------------------------------
      RETURN
      END
     
