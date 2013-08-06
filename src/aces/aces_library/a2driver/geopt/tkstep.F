      SUBROUTINE TKSTEP(GRD, HES, SCRATCH, HESMOD, GRDMOD, STEP, 
     &                  GRDHES, DIAGHES, VEC, RFAMAT, EVRFAMAT,
     &                  BMATRIX, HES_INTACT, FSCR, WORK, PMAT, 
     &                  IADITNL)
C
C Control all the optimiztion alogrithms. The old EFOL (John Stanton)
C subroutine was not flexible enough to add new optimization 
C algorithms. This routine act as a front end, and call appropriate
C routines depending on user's choice of a particular algorithm.
C
C Ajith Perera, March, 2000
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION LMBDAN, LMBDAP
      LOGICAL MORSE, TS, NRORMANR, RFA, EVFTS, IGTS, QSD,
     &        QSTLST_CLIMB
C
#include "mxatms.par"
#include "machsp.com"
#include "coord.com"
#include "jodaflags.com"
#include "flags.h"

      PARAMETER (STPTOL = 1.0D-12, LUOUT = 6,
     &           LUDONE = 80,
     &           EPS = 1.0D-8)

      INTEGER CONSTRAINS

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
C
      COMMON /OPTCTL/ IPRNT, INR, IVEC, IDIE, ICURVY, IMXSTP, ISTCRT, 
     &                IVIB, ICONTL, IRECAL, INTTYP, IDISFD, IGRDFD,
     &                ICNTYP,ISYM, IBASIS, XYZTol
C
      DIMENSION HES(NXM6,NXM6), GRD(NXM6), SCRATCH(NX*NX),
     &          HESMOD(NOPT,NOPT), GRDMOD(NOPT), STEP(NXM6), 
     &          IQFIX(3*MXATMS,3*MXATMS), DIAGHES(NOPT, NOPT), 
     &          GRDHES(NOPT),VEC(NOPT),  EVRFAMAT(NOPT+IADITNL,
     &          NOPT+IADITNL), RFAMAT(NOPT+IADITNL, NOPT+IADITNL),
     &          BMATRIX(NXM6*NX), HES_INTACT(NX*NX), FSCR(NXM6*NX),
     &          ENERGY(2),WORK(6*NX),PMAT(NX*NX)
C
      DIMENSION CONSTRAINS(3*MXATMS)
      
      Write(6,*) iFlags2(h_IFLAGS2_geom_opt)
C
C Form modified Hessian containing optimized coordinates
C only, and transform it to totally symmetric coordinates. Do the 
C same for gradient vector as well. Do all other initilizations
C and printing. Note that STPMAX is a user defined flag. It is
C initially set in mkvmol.F but is reprocessed in setup.f.
C
      IBREAK = 0 
      IF (iFlags2(h_IFLAGS2_geom_opt) .eq. 1) THEN
C
C This is a user defined internal optimization
C
         Write(6,*) "Entering setup"
         CALL SETUP(HES, GRD, HESMOD, GRDMOD, IQFIX, STPMAX, LUOUT)

      ELSE IF (iFlags2(h_IFLAGS2_geom_opt) .eq. 2) THEN
C
C This is a Cartesian optimization
C
         CALL XDCOPY(NOPT, GRD, 1, GRDMOD, 1)
C
C Vibrations and rotation must be projected out from the Hessian. In 
C this particular instance unit masses at the points where atoms are
C located. 
C
         Do Iatoms = 1, NX/3
            Fscr(Iatoms) = 1.0D0
         Enddo 
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "@-TKSTEP, The Hessian before the projection"
      call output(Hes, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
         CALL PROJEC_FC(Q, Hes, Fscr, GrdMOD, HesMOd, Scratch, 
     &                  EPS, Nx/3, .True., .True., .False.)
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "@-TKSTEP The Hessian after the projection"
      call output(Hes, 1, NX, 1, NX, NX, NX, 1)
#endif
C
C If the optimization is constrained lets handle it here instead
C of in projec_fc.F in libr (projec_fc.F is shared by other modules).
C Ajith Perera, 06/2013.

      IF (iFlags2(h_IFLAGS2_opt_control) .EQ. 3) Then

         CALL READ_CART_CONSTRAINS(CONSTRAINS, NUM_CONSTRAINS,
     &                             NOT_FOUND)
         IF (NOT_FOUND) THEN
            WRITE(6, "(1x, 2a)") "Constrained optimization requested "
     &                         , "but no constrains found. Proceeded "
            WRITE(6, "(1x,a)") "to do a full optimization."
         ELSE
            CALL APP_CART_CONSTRAINS(HES, SCRATCH, PMAT,
     &                               GRD, GRDMOD, CONSTRAINS,
     &                               NUM_CONSTRAINS, NX)
         ENDIF
      ENDIF
C
C Note that for Cartesian opts., setup call really do nothing,
C but it need to be called to update several variable, especially
C the ncycle. 
C 
         CALL SETUP(HES, GRD, HESMOD, GRDMOD, IQFIX, STPMAX, LUOUT)
         CALL XDCOPY(NOPT*NOPT, HES, 1, HESMOD, 1)
C
      ELSE IF (iFlags2(h_IFLAGS2_geom_opt) .gt. 2) THEN
C
C This is a redundant internal optimization
C
         CALL  PROJEC_IFC(HES, FSCR, SCRATCH, PMAT, HES_INTACT, 
     &                    DIAGHES, GRD, GRDMOD, NXM6)
         CALL SETUP(HES, GRD, HESMOD, GRDMOD, IQFIX, STPMAX, LUOUT)
      ENDIF
C
C Diagonalize the symmetrized Hessian. First, save the
C original symmetrized Hessian in HES to be used in QST/LST
C procedure. See the dependents of ANLYSHES that deal
C with the QST/LST procedure. Note that HESMOD and DIAGHES

C
      IF (IFLAGS2(h_IFLAGS2_opt_control) .GT. 0) CALL 
     &    XDCOPY(NOPT*NOPT, HESMOD, 1, HES, 1)

#ifdef _DEBUG_LVL0
      If (IFLAGS2(h_IFLAGS2_opt_control) .GT. 0) then
      Write(6,*) "@-TKSTEP The Hessian restored for QST/LST"
      CALL OUTPUT(HES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
      endif 
#endif
#ifdef _DEBUG_LVLM1
      Write(6,*) NX, NOPT
      Write(6,*) "@-TKSTEP The Hessian before diagonalization"
      call output(HesMod, 1, NX, 1, NX, NX, NX, 1)
#endif
      CALL EIG(HESMOD, DIAGHES, NOPT, NOPT, 1)
      IF (iFlags2(h_IFLAGS2_geom_opt) .eq. 2) CALL EVEC_SHIFT(HESMOD,
     &           DIAGHES, HES_INTACT, NOPT)
C
C Analyse the eigenvalues and eigenvectors of the current Hessian
C to decide the direction to follow. Note that ANALYSHES does
C meaningful work only for the transition state searches.
C
C ANLYSHES was modified to incorporate linear and quadratic
C synchronous (LST & QST) methods(?) of Schlegel et al. (Israel
C Journal of Chemistry, vol. 33, pg. 449, 1993). Both LST & QST
C approaches control the direction of the steps being taken until
C the quadratic region of the transition state is reached. Then
C the step direction is along the vector that corresponds to the
C smallest eigenvalue of the Hessian or a user-specified
C eigenvector of the Hessian. See ANALYSHES for further comments.
C Ajith Perera, 07/04.
C
C Important Note: RFAMAT and EVRFAMAT are used in ANLYSHES to
C                 keep the LST and QST directional vectors. They
C                 revert back to zero after their use in ANLYSHES so that
C                 there intended use in GETLAMBDA is not upset.
C
      Write(6,*) "Enetering ANLYSHES"
      CALL ANLYSHES(HESMOD, DIAGHES, RFAMAT, EVRFAMAT, HES,
     &              GRDMOD, SCRATCH, PMAT, VEC, NOPT, NX,
     &              NCYCLE, INR, IVEC, IMODE, IDIE, TS,
     &              NRORMANR, RFA, EVFTS, IGTS, QSD, LUOUT,
     &              QSTLST_CLIMB, NATOMS, IPRNT)
C
C
C Calculate the gradients along Hessian eigenvectors (GRDHES). This
C is all you need to do bare NR (MANR) optimizations. Also fill the
C scratch vector with the internal coordinates of the previous step
C (first cycle it is identical to the user given starting values!)
C
      CALL ZERO(GRDHES, NOPT)
      DO 1180 I = 1, NOPT
C
         SCRATCH(I) = R(NOPTI(I))
C
         DO 1179 J = 1, NOPT
            GRDHES(I)=GRDHES(I)+GRDMOD(J)*DIAGHES(J,I)
 1179    CONTINUE
 1180 CONTINUE
#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,*) "The gradients along eiegenvectors"
        Write(6,"(3F13.8)") (GRDHES(I), I=1, NOPT)
        Write(6,*)
#endif
C
C Do the Newton-Raphson or Morse adjusted Newton-Raphson minima search
C (Only applies to internal coordinate optimizations, Ajith Perera, 12/08).
C     
      MORSE = (INR .EQ. 3 .OR. INR .EQ. 5) .and. 
     &        (iFlags2(h_IFLAGS2_geom_opt) .ne. 2)
C
CSSS      Write(6,*) "Debug Info: The QST and LST climb", QSTLST_CLIMB
C
#ifdef _NOSKIP
      IF (NCYCLE .GE. 2 .AND. .NOT. TS) THEN
C
C The following stuff is for the line search. God I need some help
C from  litrature in order to finish this. Let's comment this out
C until, I have some help from others.
C
          CALL DGETREC(20, 'JOBARC', 'TOTENERG', 1, CURR_ENRG)
          CALL DGETREC(20, 'JOBARC', 'OLDENERG', 1, PREV_ENRG)
          CALL DGETREC(20, 'JOBARC', 'OLDGRADS', NXM6, 
     &                WORK(3*NXM6+1))
          CALL DGETREC(20, 'JOBARC', 'OLDGEOMT', NXM6,
     &                WORK(NXM6 +1))
C
          CALL XDCOPY(NXM6, GRDMOD, 1, WORK(2*NXM6+1), 1)
          CALL XDCOPY(NXM6, R,      1, WORK(1),        1)
C
          ENERGY(1) = PREV_ENRG
          ENERGY(2) = CURR_ENRG 
C
#ifdef _DEBUG_LVL0
          write(6,*)
          Write(6,*) "The previous gradient"
          Write(6,"(3F10.5)") (WORK(3*NXM6+i), i=1, NOPT)
          Write(6,*) "The current gradient"
          Write(6,"(3F10.5)") (WORK(2*NXM6+i), i=1, NOPT)
          Write(6,*) "The previous step increment"
          Write(6,"(3F10.5)") (WORK(NXM6+i), i=1, NOPT)
          Write(6,"(1x,a,2F10.5)") "Previous and current energies", 
     &                          PREV_ENRG, CURR_ENRG
#endif  
          CALL LN_SEARCH(ENERGY, WORK(1), WORK(NXM6+1), 
     &                   WORK(2*NXM6+1), WORK(3*NXM6+1),
     &                   WORK(4*NXM6+1), NXM6, EPS)
C
CC Please zero out the SCRATCH before leave this block so what's done
CC here can not interfere with the rest. 
C
      ENDIF
#endif
C
      IF (NRORMANR) THEN         
C
         CALL NEWRAPH(SCRATCH, GRDHES, HESMOD, DIAGHES, MORSE, NOPT,  
     &                NX, LUOUT)
C         
      ELSE IF (TS .OR. RFA) THEN

#ifdef _DEBUG_LVL0
      Write(6,"(a,l)") "Debug Info: The QST and LST climb",QSTLST_CLIMB
#endif
C
C First get the Lambda(P) and Lambda(N) parameters: Jon Baker
C J. Comput. Chem. 7, 385, 1986 and references their in. In general
C this is an eigenvector following method based on Rational Function Approximation
C that can apply for both minimas and transition states. In the case of minima 
C search Lambda(N) correspond to the lowest eigenvalue of RFA matrix. There
C is no Lambda(p) parameter for minima search.
C
         CALL GETLAMBDA(HESMOD, GRDHES, RFAMAT, NOPT, EVRFAMAT, 
     &                  LMBDAN, LMBDAP, IADITNL, IMODE, LUOUT, 
     &                  TS, RFA, QSTLST_CLIMB)
C         
         IF (RFA) CALL FLOWEVFMIN(SCRATCH, GRDHES, HESMOD, DIAGHES, 
     &                           LMBDAN, MORSE, NX, NOPT)
C
         IF (EVFTS) CALL FLOWEVFTS(SCRATCH, GRDHES, HESMOD, DIAGHES,
     &                             LMBDAN, LMBDAP, STPMAX, MORSE, 
     &                             IMODE, NX, LUOUT, IBREAK, NOPT,
     &                             QSTLST_CLIMB)
C         
      ELSE IF (IGTS .OR. QSD) THEN
C     
         CALL QSDMINORTS(SCRATCH, GRDMOD, HESMOD, DIAGHES, RFAMAT, 
     &                 EVRFAMAT, STPMAX, EPS, IGTS, QSD, NX, NOPT,
     &                 NCYCLE)
C
      ENDIF
C    
C Process the step taken: Convert to internal from symmetry coodinates,
C filter, apply the requested step size control, update the Geometry
C vector (R).
C
C Initialize the starting trust radius. At this point it is tied to
C the maximum step size, controlled by the user by setting MAX_STEP.
C STPMAX is set in mkvmol.F and setup.F.
C
      BEGN_TRUST_RAD = STPMAX
      CALL PROCESTEP(SCRATCH, HES, GRDMOD, SCALE, STPMAG, STEP,
     &               STPTOL, STPMAX, IQFIX, ISTCRT, BEGN_TRUST_RAD,
     &               EPS, QSTLST_CLIMB, TS)
C
C Let's do the printing. Test to see if calculation has converged.
C
      CALL SUMMARY(SCRATCH, RFAMAT, GRD, SCALE, STPMAG, IQFIX,
     &             NOPT, NX, NXM6, IBREAK, ICONTL, LUOUT, NCYCLE,
     &             LUDONE, BMATRIX, HES_INTACT, FSCR, VEC,
     &             STEP)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The R vector before leaving tkstep"
      Write(6,"(5F13.7)"), (R(I), I =1, NXM6)
#endif

C
      RETURN
      END
