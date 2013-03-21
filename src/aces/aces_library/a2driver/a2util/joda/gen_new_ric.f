      SUBROUTINE GEN_NEW_RIC(CARTCOORD, REDUNCO, IATOMICNMBER, 
     &                       NRATMS, TOTNOFBND, NW_TOTREDNCO, 
     &                       TOTNOFANG, TEST_4CHANGE)
C
C Generate the redundent internal coordinates. 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


C 
      PARAMETER(THRESHOLD = 5.0D00, EPSILON = 1.0D-10)
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH,
     &        FRAGSCR
      LOGICAL TEST_4CHANGE
C
C The following arrays need to be managed dynamically:
C
C BNDLENGTHS: Keep the bond lengths of atom pairs. It is
C             of length MXATMS*MXATMS
C IATOMICNMBER: Keep the atomic numbers of all atoms and
C               it is MAXATMS long.
C SMOFCOVRAD: The sum of covalent radius of pairs of atoms
C             and of MXATMS*MXATMS long.
C NCONPRCNTR: The number of connectivities for an atom. The
C             maximum number is set to 10, so this is a
C             rather small (it is not currently being
C             used, but needs in the future).
C
      DIMENSION CARTCOORD(3*MXATMS), BNDLENTHS(MXATMS, MXATMS), 
     &          IATOMICNMBER(MXATMS), SMOFCOVRADI(MXATMS, MXATMS),
     &          IBNDTO(MXATMS*MXATMS), NCONPRCNTR(MXATMS),
     &          IREDUNCO(4, MAXREDUNCO), 
     &          REDUNCO(MAXREDUNCO), MARK_FRAGMENTS(MXATMS, MXATMS),       
     &          LENGTH_FRAGEMENTS(MXATMS), VECBA(3), VECBC(3), 
     &          VECBAD(3), VECCB(3), VECAB(3), VECCD(3), VECABC(3),
     &          VECBCD(3),
     &          FRAGSCR(MXATMS*MXATMS),BLNGT_IJFRGS(MXATMS)
C      
      DATA  IZERO, MONE, ZERO /0, -1, 0.0D0/
C
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0 
      PI     = (ATAN(DFLOAT(1))*DFLOAT(4))
C
      NW_TOTREDNCO = IZERO
      BHORTOANG    = 0.529177249D0
C
      CALL GETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
      CALL GETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO,
     &            IREDUNCO)   
      CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
      CALL GETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
      CALL GETREC(20, 'JOBARC', 'IBONDTO ', NRATMS*NRATMS,
     &            IBNDTO)
C
      Print*, "The redundent internal coordinates assignments"
      Write(6,*)
      Do i = 1, TOTREDNCO
         Write(6,111) (iredunco(j, i), j=1, 4)
      Enddo
C
C Assign bond coordiante
C

      DO IBNDS = 1, TOTNOFBND
         IF (IREDUNCO(2, IBNDS) .NE. 0) THEN
            ICON1  = IREDUNCO(1, IBNDS)
            ICON2  = IREDUNCO(2, IBNDS)
            DISTAB = DIST(CARTCOORD(3*ICON1 - 2), 
     &               CARTCOORD(3*ICON2 - 2))
            REDUNCO(IBNDS) = DISTAB
C     
            WRITE(6,"(a,F10.5)") "The bond distance =",DISTAB*0.529177249d0
         ENDIF
      ENDDO
C
C Assign bond angle Coordinates.
C
      DO IANGS = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
     
         IF (IREDUNCO(4, IANGS) .EQ. MONE) THEN
C
C This is for linear arrangements (two angle coordinates are
C required). 
C
            REDUNCO(IANGS) = 180.0D0*DINVPI
            WRITE(6,"(a,F10.5)") "The Bond Angle =", ANGL/DINVPI
C
         ELSE
C
C This is for the non-linear arrangements.
C
            ICON1  = IREDUNCO(1, IANGS)
            ICON2  = IREDUNCO(2, IANGS)  
            ICON3  = IREDUNCO(3, IANGS)
            CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2),
     &               VECBC, 1)
            CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON1 - 2),
     &               VECBA, 1)
            ANGL    = ANGLE(VECBC, VECBA, 3)*DINVPI
            REDUNCO(IANGS) = ANGL
            WRITE(6,"(a,F10.5)") "The Bond Angle =", ANGL/DINVPI
         ENDIF
C
      ENDDO
C
C Assign bond dihedral angle Coordinates.
C
      DO IDIHS = (TOTNOFANG + TOTNOFBND + 1),  TOTREDNCO
         ICON1  = IREDUNCO(1, IDIHS)
         ICON2  = IREDUNCO(2, IDIHS)
         ICON3  = IREDUNCO(3, IDIHS)
         ICON4  = IREDUNCO(4, IDIHS)
         CALL VEC(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2),
     &         VECAB, 1)
         CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2),
     &         VECBC, 1)
         CALL VEC(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2),
     &         VECCB, 1)
         CALL VEC(CARTCOORD(3*ICON4 - 2), CARTCOORD(3*ICON3 - 2),
     &         VECCD, 1)
C
         DISTAB = DIST(CARTCOORD(3*ICON1 - 2), 
     &            CARTCOORD(3*ICON2 - 2))
         DISTBC = DIST(CARTCOORD(3*ICON3 - 2), 
     &            CARTCOORD(3*ICON2 - 2))
         DISTCD = DIST(CARTCOORD(3*ICON3 - 2), 
     &            CARTCOORD(3*ICON4 - 2))
C
C First evaluate the dihedral angle. This is calculated as the
C angle between the two unit vectors that are perpendicular to
C the ABC and BCD planes for A-B-C-D pattern.
C
         CALL CROSS(VECAB, VECBC, VECABC, 1)
C
C We need the vector CD, obtain that from vec DC (note the misnomar)
C
         CALL DSCAL(3, -1.0D0, VECCD, 1)
         CALL CROSS(VECBC, VECCD, VECBCD, 1)
C
         DANG = (ANGLE(VECABC, VECBCD, 3))*DINVPI
C
C--- Obtainig the sign of the dihedral angle:
C
         CALL DSCAL(3, -1.0D0, VECCD, 1)
         CALL CROSS(VECCD, VECCB, VECBCD, 0)
C
         SENSE = DDOT(3, VECAB, 1, VECBCD, 1)

         IF (SENSE .GT. 0.0D0) DANG = -DANG
         REDUNCO(IDIHS) = DANG
         
C
C
      ENDDO 
C
      IF (.NOT. TEST_4CHANGE) THEN
          NW_TOTREDNCO = TOTREDNCO
          RETURN
      ENDIF
C
C Assign new set of redundent internal coordiantes to check whether
C the number of redundent internal coordiantes are changing.
C
      CALL ASIGN_CNTVTES(CARTCOORD, BNDLENTHS, IATOMICNMBER,
     &                   SMOFCOVRADI, NCONPRCNTR, IBNDTO, MAXCNTVS,
     &                   NRATMS, MARK_FRAGMENTS, LENGTH_FRAGEMENTS,
     &                   FRAGSCR, BLNGT_IJFRGS)
      Print*, "Entering assign bonds"
      Write(6,*)
      CALL ASSIGN_BONDS(IBNDTO, IREDUNCO, NW_TOTREDNCO, TOTNOFBND,
     &                  NRATMS, MAXREDUNCO)
      Print*, "Entering assign angles"
      Print*, "Total # of bonds:", TOTNOFBND
      CALL ASSIGN_ANGLS(CARTCOORD, IBNDTO, IREDUNCO, NW_TOTREDNCO,
     &                  TOTNOFBND, TOTNOFANG, NRATMS, MAXREDUNCO,
     &                  THRESHOLD)
      Print*, "Entering assign dihedral angles"
      Print*, "Total # of Angles:", TOTNOFANG
      CALL ASSIGN_DIHLS(IBNDTO, IREDUNCO, NW_TOTREDNCO, TOTNOFBND,
     &                  TOTNOFANG, NRATMS, TOTNOFDIH, MAXREDUNCO)

      Write(6,*)
      Print*, "Total # of dihedrals:", TOTNOFDIH
      Print*, "The total redundent internal", TOTREDNCO
      Write(6,*)
      Print*, "The redundent internal coordinates assignments"
      Write(6,*)
      Do i = 1, TOTREDNCO
         Write(6,111) (iredunco(j, i), j=1, 4)
      Enddo
 111  Format(5X, 4(I3, 1X))
      IF (NW_TOTREDNCO .NE. TOTREDNCO) THEN
          Write(6,*) 
          Write(6,"(1x,a,a,a,a)")"The total number of redundent ",
     &                           "internals have changed during ", 
     &                           "the optimization, please start ",
     &                           "from the following coords (Angs.)"
C
          IOFF = 1
          DO I =1, NRATMS
             WRITE(6,100)ZSYM(I)(1:2),(BHORTOANG*CARTCOORD(J),
     &                                 J=IOFF,IOFF+2)
             IOFF=IOFF+3
          END DO


         CALL ERREX
      ENDIF 
 100  FORMAT(T6,A,T10,F14.8,T25,F14.8,T40,F14.8)

      RETURN
      END

