
      SUBROUTINE EVAL_QSTLST_PATH_RIC(QST_TANGENT, LST_TANGENT, WORK,
     &                                PMAT, LST, QST)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

#include "mxatms.par"
#include "machsp.com"
#include "coord.com"

      LOGICAL XYZIN, NWFINDIF, QST, LST
      DOUBLE PRECISION LST_NORM_SQR, LST_NORM_INV, LST_TANGENT
      INTEGER TOTREDNCO

      COMMON /USINT/NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/XYZIN,NWFINDIF

      DIMENSION WORK(NX*NX), LST_TANGENT(*), QST_TANGENT(*),
     &          PMAT(NX*NX)
C
C Use Eqn. 1 and 8 in Schlegel to evaluate the QST and LST normalized
C tangent vectors.
C
C Read the guess structures (user provided) for the transition
C state and the products.
C 
      CALL IGETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)

      IRECTNT = 1
      IPRDUCT = IRECTNT + TOTREDNCO
      ICURENT = IPRDUCT + TOTREDNCO
      ISCR    = ICURENT + TOTREDNCO
      IEND    = ISCR    + TOTREDNCO

      IF (IEND .GT. NX*NX) CALL INSMEM("@EVAL_QSTLST_PATH_RIC",
     &                                  IEND, NX*NX)

C Cautionary note: The reactant coordinates are copied to the
C memory location desiganted as IPRDUCT. This was as a result
C of a earlier confusion in designating multiple structures. This
C has no effect on actual work done, Ajith Perera, 12/2012

      CALL DGETREC(20,'JOBARC','RX_RICS_',TOTREDNCO,
     &            WORK(IPRDUCT))

CSSS      CALL CNVRT2RAD(WORK(IPRDUCT),3*NATOMS)
CSSS      CALL SQUSH(WORK(IPRDUCT),3*NATOMS)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Debug-Info: The Reactant Structure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I),I= 0,TOTREDNCO-1)
#endif

      IF (LST) THEN
        CALL DGETREC(1,'JOBARC','PR_RICS_',TOTREDNCO,
     &              WORK(IRECTNT))

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Debug-Info: Product Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)
#endif
        CALL XDAXPY(TOTREDNCO, -1.0D0, WORK(IRECTNT), 1,
     &             WORK(IPRDUCT), 1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Debug-Info: The Diff. Structure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I), I= 0,TOTREDNCO-1)
#endif

CSSS            CALL XDCOPY(NX, WORK(IPRDUCT), 1, LST_TANGENT, 1)
CSSS            CALL GETREC(20, "JOBARC", "BMATRIX ", 3*NATOMS*NOPT*
CSSS     &                  IINTFP, WORK(1))
CSSS            IOFFSET = 3*NATOMS*NOPT + 1
CSSS            CALL XGEMM("N", "N", NOPT, 1, 3*NATOMS, 1.0D0, WORK,
CSSS     &                  NOPT, LST_TANGENT, 3*NATOMS, 0.0D0, WORK
CSSS     &                  (IOFFSET), NOPT)
CSSS            CALL XDCOPY(NOPT, WORK(IOFFSET), 1, LST_TANGENT, 1)
CSSS        CALL SYMUNQONLY(WORK(IPRDUCT), LST_TANGENT)

        CALL XGEMM("N", "N", TOTREDNCO, 1, TOTREDNCO, 1.0D0, 
     &              PMAT, TOTREDNCO, WORK(IPRDUCT),
     &              TOTREDNCO, 0.0D0, WORK(ISCR), TOTREDNCO)
        CALL XDCOPY(TOTREDNCO, WORK(ISCR), 1, WORK(IPRDUCT), 1)

        CALL NORMAL(WORK(IPRDUCT), TOTREDNCO)
        TNORM = XDDOT(TOTREDNCO, WORK(IPRDUCT),1,WORK(IPRDUCT), 1)
C
        IF (TNORM .NE. 1.0D0) THEN
           WRITE(6,*)
           WRITE(6,"(a,a)") "Norm condition for the LST tangent is ",
     &                      "incorrect! Can not proceed."
            CALL ERREX
        ENDIF
C
        CALL XDCOPY(TOTREDNCO, WORK(IPRDUCT), 1, LST_TANGENT, 1)

#ifdef _DEBUG_LVL0
            Write(6,"(a)") "Debug-Info: LST tangent"
            Write(6,"(3(1x,F12.6))") (LST_TANGENT(I), I= 1, NOPT)
#endif

      ELSE IF (QST) THEN

C Cautionary note: The reactant coordinates are copied to the
C memory location desiganted as IRECTNT. This was as a result
C of a earlier confusion in designating multiple structures. This
C has no effect on actual work done, Ajith Perera, 12/2012

        CALL DGETREC(1,'JOBARC','PR_RICS_',TOTREDNCO,
     &              WORK(IRECTNT))

CSSS        CALL CNVRT2RAD(WORK(IRECTNT), 3*NATOMS)
CSSS        CALL SQUSH(WORK(IRECTNT), 3*NATOMS)
C
        CALL XDCOPY(TOTREDNCO, R, 1, WORK(ICURENT), 1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Debug-Info: Product Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)
      Write(6,"(a)") "Debug-Info: Reactant Structure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I), I= 0,TOTREDNCO-1)
      Write(6,"(a)") "Debug-Info: Current Structure"
      Write(6,"(3(1x,F12.6))") (Work(ICURENT + I), I= 0,TOTREDNCO-1)
#endif

        CALL XDAXPY(TOTREDNCO, -1.0D0, WORK(ICURENT), 1,
     &             WORK(IPRDUCT), 1)
        CALL XDAXPY(TOTREDNCO, -1.0D0, WORK(ICURENT), 1,
     &             WORK(IRECTNT), 1)

        QST_NORM_SQR_PX = XDDOT(TOTREDNCO, WORK(IPRDUCT), 1,
     &                         WORK(IPRDUCT), 1)
        QST_NORM_SQR_RX = XDDOT(TOTREDNCO, WORK(IRECTNT), 1,
     &                         WORK(IRECTNT), 1)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Debug-Info: P-X tructure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I), I= 0,TOTREDNCO-1)
      Write(6,"(a)") "Debug-Info: R-X Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)
#endif
        QST_NORM_SQR_PX = XDDOT(TOTREDNCO, WORK(IPRDUCT), 1,
     &                         WORK(IPRDUCT), 1)
        QST_NORM_SQR_RX = XDDOT(TOTREDNCO, WORK(IRECTNT), 1,
     &                         WORK(IRECTNT), 1)
        QST_NORM_PX_RX  = XDDOT(TOTREDNCO, WORK(IRECTNT), 1,
     &                         WORK(IPRDUCT), 1)

        QST_NORM_PX = DSQRT(QST_NORM_SQR_PX)
        QST_NORM_RX = DSQRT(QST_NORM_SQR_RX)

        IF ((QST_NORM_SQR_PX.LE.1.0D-14).OR.
     &     (QST_NORM_SQR_RX.LE.1.0D-14)    ) THEN
          Write(6,*)
          Write(6,"(a)")"AN INPUT ERROR: INSPECT THE INPUT STRUCTURES"
          CALL ERREX
        END IF

        QST_NORM_PX_INV = 1.0D0/(QST_NORM_SQR_PX)
        QST_NORM_RX_INV = 1.0D0/(QST_NORM_SQR_RX)

        CALL XDSCAL(TOTREDNCO, QST_NORM_PX_INV, WORK(IPRDUCT), 1)
        CALL XDSCAL(TOTREDNCO, QST_NORM_RX_INV, WORK(IRECTNT), 1)

        CALL XDAXPY(TOTREDNCO, -1.0D0, WORK(IRECTNT), 1,
     &             WORK(IPRDUCT), 1)
        DENOM = (QST_NORM_SQR_PX + QST_NORM_SQR_RX - 
     &           2.0D0*QST_NORM_PX_RX)

        IF (DENOM .LE. 1.0D-14) THEN
           Write(6,*)
           Write(6,"(a)")"AN INPUT ERROR: INSPECT THE INPUT STRUCTURES"
           CALL ERREX
        ENDIF

        A_SCALE = DSQRT((QST_NORM_SQR_PX*QST_NORM_SQR_RX)/DENOM)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Various Norms"
      Write(6,"(a,2(1x,F12.6))") "QST_NORM_SQR_PX,RX :",QST_NORM_SQR_PX,
     &                          QST_NORM_SQR_RX
      Write(6,"(a,2(1x,F12.6))") "QST_NORM_PX,RX     :",QST_NORM_PX,
     &                          QST_NORM_RX
      Write(6,"(a,(1x,F12.6))") "QST_NORM_PX_RX     :",QST_NORM_PX_RX
#endif
CSSS            CALL XDCOPY(3*NATOMS, WORK(IPRDUCT), 1, QST_TANGENT, 1)
CSSS            CALL GETREC(20, "JOBARC", "BMATRIX ", 3*NOPT*NATOMS*
CSSS     &                  IINTFP, WORK(1))
CSSS            IOFFSET = 3*NATOMS*NOPT + 1
CSSS            CALL XGEMM("N", "N", NOPT, 1, 3*NATOMS, 1.0D0, WORK,
CSSS     &                  NOPT, QST_TANGENT, 3*NATOMS, 0.0D0, WORK
CSSS     &                  (IOFFSET), NOPT)

        CALL XDSCAL(TOTREDNCO, A_SCALE, WORK(IPRDUCT), 1)
C
        CALL XGEMM("N", "N", TOTREDNCO, 1, TOTREDNCO, 1.0D0,
     &              PMAT, TOTREDNCO, WORK(IPRDUCT),
     &              TOTREDNCO, 0.0D0, WORK(ISCR), TOTREDNCO)
        CALL XDCOPY(TOTREDNCO, WORK(ISCR), 1, WORK(IPRDUCT), 1)

        CALL NORMAL(WORK(IPRDUCT), TOTREDNCO)

        TNORM = XDDOT(TOTREDNCO, WORK(IPRDUCT),1,WORK(IPRDUCT), 1)
C
        IF (TNORM .NE. 1.0D0) THEN
           WRITE(6,*)
           WRITE(6,"(a,a)") "Norm condition for the QST tangent is ",
     &                      "incorrect! Can not proceed."
           CALL ERREX
        ENDIF
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Checking the Norm of QST tangent"
      Write(6,"(a,1x,F12.6)")"Tangent: ", TNORM
#endif
C
CSSS        CALL SYMUNQONLY(WORK(IPRDUCT), QST_TANGENT)
          
         CALL XDCOPY(TOTREDNCO, WORK(IPRDUCT), 1, QST_TANGENT, 1)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Debug-Info: QST tangent"
      Write(6,"(3(1x,F12.6))") (QST_TANGENT(I), I= 1, NOPT)
#endif
         END IF

       RETURN
       END

