      SUBROUTINE CONVQ(A,FC,FI,IX,ISYM)
C
C CONVERTS FORCES IN CARTESIAN COORDINATES TO INTERNAL COORDINATES (IX=0
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL XYZIN, NWFINDIF
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)









































































































































































































































































































































































































































































































































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



      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
      DIMENSION A(3*NATOMS,NXM6),FC(NX),FI(NXM6)
C
C
c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
C
      IF(IX.EQ.0)THEN
         DO 10 I=1,NXM6
            Z=0.D0
            DO 20 J=1,NX
               Z=FC(J)*A(J,I)+Z
 20         CONTINUE
            FI(I)=Z
 10      CONTINUE

C

         Write(6,*)
         Write(6,*) "@-CONVQF,  The outgoing internal gradients"
         WRITE(6,"(3F13.7)") (FI(IJL), IJL = 1, NXM6)
         Write(6,*)
C
C
C SYMMETRIZE THE GRADIENT IF NECESSARY
C
         IF (ISYM.EQ.1)THEN
            IF (XYZIN) THEN
               IF (IFLAGS(60).EQ. 2) THEN
                  CALL GETREC(20, 'JOBARC', 'PLSMINSP', NXM6, NCON)
               END IF
               DO I = 1, NUNIQUE
                  IP = IUNIQUE(I)
cSSS                  WRITE(6,*) "The IP =", IP, NEQ(IP),FI(IP) 
                  Z  = FI(IP)
C
                  DO J = 1, NEQ(IP)
                     IF (NCON(IEQUIV(I, J)) .EQ. 1) THEN
                        SIGN = -1.0D0
                     ELSE
                        SIGN =  1.0D0
                     ENDIF
cSSS                     WRITE(6,*) FI(IEQUIV(I,J)), IEQUIV(I,J),SIGN
                     Z = Z + SIGN*FI(IEQUIV(I, J))
                  ENDDO
C
                  FIAVE=Z/(NEQ(IP)+1)
cSSS                  WRITE(6,*) FIAVE
                  DIFF=DABS(FI(IP)-FIAVE)
                  IF(DIFF.GT.1D-7)WRITE(LuOut,72)DIFF,IP
C
 72               FORMAT(T3,'@CONVQ-W, Nominally equivalent',
     $            'internal coordinate gradients vary by at least ',
     $            e10.5,'.'/
     $            T3,' Problem with gradients equivalent to ',
     $            'parameter [',I3,'].')
                  FI(IP) = FIAVE
C
                  DO J = 1, NEQ(IP)
                     IF (NCON(IEQUIV(I, J)) .EQ. 1) THEN
                        SIGN = -1.0D0
                     ELSE
                        SIGN = 1.0D0
                     ENDIF
                     FI(IEQUIV(I, J))= SIGN*FIAVE
cSSS                     WRITE(6,*) "The IP =", IEQUIV(I,J),FI(IEQUIV(I, J))
                  ENDDO
               ENDDO
               CALL FILTER(FI, NXM6, 1.0D-6) 
               RETURN
            ELSE
               DO 100 I=1,NUNIQUE
                  IP=IUNIQUE(I)
                  Z=FI(IP)
                  DO 102 J=1,NEQ(IP)
 102                 Z=Z+FI(IEQUIV(I,J))
                     FIAVE=Z/(NEQ(IP)+1)
                     DIFF=DABS(FI(IP)-FIAVE)
                     IF(DIFF.GT.1D-7)WRITE(LuOut,71)DIFF,IP
 71                 FORMAT(T3,'@CONVQ-W, Nominally equivalent',
     $             'internal coordinate gradients vary by at least ',
     $              e10.5,'.'/
     $              T3,' Problem with gradients equivalent to ',
     $              'parameter [',I3,'].')
                    FI(IP)=FIAVE
C
                    DO J = 1, NEQ(IP)
                       FI(IEQUIV(I,J))= FIAVE
                    ENDDO
C
 100             CONTINUE
            ENDIF
            RETURN
         ENDIF
      ELSE
         DO 30 I=1,NX
            Z=0.D0
            DO 40 J=1,NXM6
 40            Z=A(I,J)*FI(J)+Z
            FC(I)=Z
 30      CONTINUE
         RETURN
      ENDIF
      RETURN
      END
