C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      SUBROUTINE GEN_NEW_RIC(CARTCOORD, REDUNCO, IATOMICNMBER, 
     &                       NRATMS, TOTNOFBND, NW_TOTREDNCO)
C
C Generate the redundent internal coordinates. 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "io_units.par"
C 
      PARAMETER(THRESHOLD = 0.5D00, EPSILON = 1.0D-10)
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
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
      DIMENSION CARTCOORD(3*NRATMS), BNDLENTHS(NRATMS, NRATMS), 
     &          IATOMICNMBER(NRATMS), SMOFCOVRADI(NRATMS, NRATMS),
     &          IBNDTO(NRATMS*NRATMS), NCONPRCNTR(MAXCNTVS),
     &          IREDUNCO(4, MAXREDUNCO), 
     &          REDUNCO(MAXREDUNCO), MARK_FRAGMENTS(NRATMS, NRATMS),       
     &          LENGTH_FRAGEMENTS(NRATMS), VECBA(3), VECBC(3), 
     &          VECBAD(3), VECCB(3), VECAB(3), VECCD(3), VECABC(3),
     &          VECBCD(3)
C      
      DATA  IZERO, MONE, ZERO /0, -1, 0.0D0/
C
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0 
      NW_TOTREDNCO = IZERO
C
#ifdef _NO_SKIP
C
C Reassigning connectivites during iterative determination of 
C Cartesians from internals can be problematic, especially if
C isolated fragments are invloved. If someone wants to be 
C adventurous, define _NOSKIP during compilation.

#ifdef _DEBUG_LVLM1
      Print*, "Entering assign connectivities"
      Print*, (IATOMICNMBER(I), I=1, NRATMS)
#endif
C
      CALL ASIGN_CNTVTES(CARTCOORD, BNDLENTHS, IATOMICNMBER, 
     &                   SMOFCOVRADI, NCONPRCNTR, IBNDTO, MAXCNTVS,
     &                   NRATMS, MARK_FRAGMENTS, LENGTH_FRAGEMENTS)
C
#else
C
      CALL IGETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
      CALL IGETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO,
     &            IREDUNCO)   
      CALL IGETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL IGETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
      CALL IGETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
      CALL IGETREC(20, 'JOBARC', 'IBONDTO ', NRATMS*NRATMS,
     &            IBNDTO)
#endif
C
#ifdef _DEBUG_LVL0 
      Print*, "The redundent internal coordinates assignments"
      Write(6,*)
      Do i = 1, TOTREDNCO
         Write(6,111) (iredunco(j, i), j=1, 4)
      Enddo
 111  Format(5X, 4(I3, 1X))
      Write(6,*) 
      Write(6,"(a)") "The Ibndto array"
      Write(6,"(8I4)") (Ibndto(I), I=1, nratms*nratms)
      Write(6,*)
#endif
C 
C Assign bond length Coordinates.
C
#ifdef _NO_SKIP
      CALL ASSIGN_BONDS(IBNDTO, IREDUNCO, NW_TOTREDNCO, TOTNOFBND,
     &                  NRATMS, MAXREDUNCO)
C
#endif
      DO IBNDS = 1, TOTNOFBND
         IF (IREDUNCO(2, IBNDS) .NE. 0) THEN
            ICON1  = IREDUNCO(1, IBNDS)
            ICON2  = IREDUNCO(2, IBNDS)
            DISTAB = DIST(CARTCOORD(3*ICON1 - 2), 
     &               CARTCOORD(3*ICON2 - 2))
            REDUNCO(IBNDS) = DISTAB
         ENDIF
      ENDDO
C
C Assign bond angle Coordinates.
C
#ifdef _NO_SKIP
      CALL ASSIGN_ANGLS(CARTCOORD, IBNDTO, IREDUNCO, NW_TOTREDNCO,
     &                  TOTNOFBND, TOTNOFANG, NRATMS, MAXREDUNCO,
     &                  THRESHOLD)
#endif
C
      DO IANGS = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
     
         IF (IREDUNCO(4, IANGS) .EQ. MONE) THEN
C
C This is for linear arrangements (two angle coordinates are
C required). 
C
            REDUNCO(IANGS) = ZERO
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
         ENDIF
C
      ENDDO
C
C Assign bond dihedral angle Coordinates.
C
#ifdef _DEBUG_LVL0
      Print*, "Entering assign dihedral angles"
      Print*, "Total # of Angles:", TOTNOFANG 
#endif 
C
#ifdef _NO_SKIP
      CALL ASSIGN_DIHLS(IBNDTO, IREDUNCO, NW_TOTREDNCO, TOTNOFBND,
     &                  TOTNOFANG, NRATMS, TOTNOFDIH, MAXREDUNCO)
C
C Assign dihedral angles for nearly planar structures.
C
      IF (TOTNOFDIH .EQ. IZERO) THEN
C
#ifdef _DEBUG_LVL0
      Print*, "Entering assign planar dihedrals"
#endif 

         CALL ASSIGN_DIHLSNPL(CARTCOORD, IBNDTO, IREDUNCO,  
     &                        NW_TOTREDNCO, TOTNOFBND, TOTNOFANG, 
     &                        TOTNOFDIH, NRATMS, MAXREDUNCO, 
     &                        THRESHOLD)
      ENDIF
#endif
C
      DO IDIHS = (TOTNOFANG + TOTNOFBND + 1),  NW_TOTREDNCO
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
         DANG = (ANGLE(VECABC, VECBCD, 3))*DINVPI
         REDUNCO(IDIHS) = DANG
#ifdef _DEBUG_LVLM1
      Print*, "The vectors making the dihedral angle"
      Print*,  (CARTCOORD(3*ICON1 - 2 + I), I=0, 2)
      Print*,  (CARTCOORD(3*ICON2 - 2 + I), I=0, 2)
      Print*,  (CARTCOORD(3*ICON3 - 2 + I), I=0, 2)
      Print*,  (CARTCOORD(3*ICON4 - 2 + I), I=0, 2)
      Print*, "The vec. BC: ", (VECBC(I), I=1, 3)
      Print*, "The vec. AB: ", (VECAB(I), I=1, 3)
      Print*, "The vec. CB: ", (VECCB(I), I=1, 3)
      Print*, "The vec. CD: ", (VECCD(I), I=1, 3)
      Print*, "The angle:", DANG/DINVPI
#endif

      ENDDO 
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The total redundent internal", NW_TOTREDNCO
#endif

      RETURN
      END

