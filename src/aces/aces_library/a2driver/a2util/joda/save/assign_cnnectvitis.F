      SUBROUTINE ASIGN_CNTVTES(CARTCOORD, BNDLENTHS, IATOMICNMBER, 
     &                         SMOFCOVRADI, NCONPRCNTR, IBNDTO,
     &                         MAXCNTVS, NRATMS, ICONNECTO,
     &                         CONNECTED_ATOMS, MARK_FRAGMENTS,
     &                         LENGTH_FRAGMENTS, DIST_IJFRAGS)
C
C We need to decide some sort of a connectivity pattern if want to do
C Cartesian optimizations based on RIC. I prefer that the user specify the 
C connectivities, but Gaussian crowd certainly going to cry!
C So Let's use this covalent radii scheme Dr. Fau is talking about. 
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION PNTSEV
      LOGICAL FLAG, NOPClFPS, FRAGMENT
      INTEGER AATOM, BATOM, VALENCY, CONNECTED_ATOMS, IND, IDUM,
     &        JDUM
CSSS    NEIGHBORS(NRATMS, NRATMS)
C
      DIMENSION CARTCOORD(3*NRATMS), BNDLENTHS(NRATMS, NRATMS), 
     &          IATOMICNMBER(NRATMS), SMOFCOVRADI(NRATMS, NRATMS),
     &          IBNDTO(NRATMS, NRATMS), NCONPRCNTR(NRATMS),
     &          ICONNECTO(NRATMS, MAXCNTVS), VALENCY(86),
     &          CONNECTED_ATOMS(NRATMS), MARK_FRAGMENTS(NRATMS, 
     &          NRATMS), LENGTH_FRAGMENTS(NRATMS)
C      
      DIMENSION VRADII(17)
      DATA VRADII /1.20D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &             1.55D0, 1.52D0, 1.47D0, 0.0D0, 0.0D0, 0.0D0,
     &             0.0D0, 0.0D0, 1.80D0, 1.80D0, 1.75D0/
      DATA VALENCY /1, 0, 1, 2, 3, 4, 4, 2, 1, 0, 1, 2, 6, 6, 5,
     &              6, 1, 0, 1, 2, 6, 6, 6, 8, 6, 6, 6, 6, 6, 3,
     &              4, 3, 2, 1, 0, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6,
     &              6, 6, 3, 4, 3, 2, 1, 0, 1, 2, 12, 6, 6, 6, 
     &              6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     &              6, 6, 6, 6, 6, 3, 4, 3, 2, 1, 0, 1/
      DATA PNTSEV /0.71D0/, ZZERO /0.D0/, HALF/0.50D0/
      DATA IONE /1/, JZERO /0/
      DATA ZIONE /1.D0/
      DATA AMTOBOHR /1.8897265/
C
C$$$      INDEX(IDUM, JDUM) = (IDUM - 2)*(IDUM - 1)/2 + JDUM 
C
C Calculate internuclear distances; The calculated distances are 
C in Angstroms Call Zero sets up a 0 matrix for bond lenghts
C
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
      CALL ZERO(BNDLENTHS, NRATMS*NRATMS)
      CALL IZERO(IBNDTO, NRATMS*NRATMS)
C      
      INDEX1 = 1
      INDEX2 = 4
C
      DO 10 IBNDS = 2, NRATMS
         DO 20 JBNDS = 1, IBNDS - 1
            BNDLENTHS(JBNDS, IBNDS) = DIST(CARTCOORD(INDEX1),
     &                                CARTCOORD(INDEX2))
            INDEX1 = INDEX1 + 3 
 20      CONTINUE
         INDEX2 = INDEX2 + 3 
         INDEX1 = IONE
 10   CONTINUE
C
C Now calculate the scaled sum of the covalent radii for each unique
C pair of atoms
C     
cSSS CALL GETREC(20, 'JOBARC', 'ATOMCHRG', NRATMS, IATOMICNMBER)
C      
      CALL GETCOVLRADI(IATOMICNMBER, SMOFCOVRADI, NRATMS)
C      
C Now assign preliminary connectivites based on simple fact that
C length between two centers is less than or equal to 1.3 times sum of the 
C covalent radii of the centers. If they are connected, assign integer
C 1 for that pair, if not leave it as 0.
C This simple criteria is used by Helgaker and many others to assign 
C connectivities. 
C
      DO IBNDS = 2, NRATMS
         DO JBNDS = 1, (IBNDS - 1)
            IF (IBNDTO(JBNDS,IBNDS).EQ.JZERO.AND.
     &          BNDLENTHS(JBNDS,IBNDS).LE.SMOFCOVRADI(JBNDS,IBNDS)
     &         ) THEN
                IBNDTO(JBNDS,IBNDS) = IONE
#ifdef _DEBUG_LVL0
            Print*, "Bonded atoms, bond-length, sum of c-radii"
            Print*,jbnds, ibnds, BNDLENTHS(JBNDS, IBNDS)*0.529177249d0, 
     &             SMOFCOVRADI(JBNDS,IBNDS)*0.529177249d0, 
     &             SMOFCOVRADI(IBNDS,JBNDS)*0.529177249d0
#endif
            END IF
         END DO
      END DO
C
#ifdef _NOSKIP
C This is what I had when I implemented first. It turns out that the
C above criteria works better. At the moment, I am skipping this,
C but leaving this  the have the possiblity to use it if needed. 
C Turning on this do no harm at all.
C
      DO IBNDS = 2, NRATMS
cSSS         NCONPRCNTR(IBNDS) = JZERO
         DO JBNDS = 1, (IBNDS - 1)
C
C Calculate bond orders; The BNDORDR array is set to 
C 0, 1, 2, 3, 4, 5 for No bond, single, double, triple, hydrogen and
C bond dimer.
C
            BNDORDER = EXP((SMOFCOVRADI(JBNDS, IBNDS) -
     &                 BNDLENTHS(JBNDS, IBNDS))/PNTSEV)
C
            IF (BNDORDER .GE.ZIONE .AND. IBNDTO(JBNDS,IBNDS)
     &          .EQ. JZERO) THEN
C               IBNDTO(JBNDS, IBNDS) = MAX(IDINT(BNDORDER + HALF), 3)
                IBNDTO(JBNDS, IBNDS) = IONE
            ELSE
               IBNDTO(JBNDS, IBNDS) = JZERO
            ENDIF
         END DO
      END DO
C
C A trial at developing a procedure to locate the existence
C isolated fragaments.
C
      DO IATOMS = 2, NRATMS
        DO JATOMS =1, (IATOMS - 1)
           IF (IBNDTO(JATOMS, IATOMS) .EQ. JZERO) THEN
c              DO ILEVELS = 1,  MAXLEVELS
c                 DO JLEVLES = 1, MAXLEVELS

                    INO_CTVS = NCONPRCNTR(IATOMS)
                    JNO_CTVS = NCONPRCNTR(JATOMS)
                    IATMCNUM = IATOMICNMBER(IATOMS)
                    JATMCNUM = IATOMICNMBER(JATOMS)
#ifdef _DEBUG_LVL0
      Print*, "The atomic numbers of the unconnected pair",
     &         IATMCNUM, JATMCNUM 
#endif     

                     IF (INO_CTVS .LT. VALENCY(IATMCNUM) .AND.
     &                   JNO_CTVS .LT. VALENCY(JATMCNUM)) THEN 
C    
c                    DO WHILE (FRAGMENT)
                       ICONS = 0
                       FRAGMENT = .TRUE. 
                       DO WHILE (FRAGMENT.AND.ICONS.LT.INO_CTVS)
                          ICONS = ICONS + 1
                          AATOM = ICONNECTO(IATOMS, ICONS)
                          JCONS = 0
                          DO WHILE (FRAGMENT.AND.JCONS.LT.JNO_CTVS)
                             JCONS = JCONS + 1
                             BATOM = ICONNECTO(JATOMS, JCONS)
#ifdef _DEBUG_LVL0
      Print*, "Atoms bonded to to disconnected pair"
      Print*, "The current pair:", IATOMS, JATOMS
      Print*, "The current pair connected to:", AATOM, BATOM
#endif
                             IF (IBNDTO(AATOM, BATOM) .EQ. 1)
     &                           FRAGMENT = .FALSE.
                          END DO
                       END DO
                      ENDIF 
           Do ICTVS = 1, INO_CTVS
              AATOM = ICONNECTO(IATOMS, ICONS)
              DO LOC_AATOM =1, NATOSM`
              NLIST(LOC_AATOM) = AATOM
            

c          ! When we get here, one of two possibilities exist:
c          ! FRAGMENT=FALSE and we quit the outter-most loop OR
c          ! ICONS=INO_CTVS and JCONS=JNO_CTVS and FRAGMENT=TRUE
c
c          ! This code sample looks incomplete because nothing
c          ! changes the state of FRAGMENT and the outter-most
c          ! loop will run forever.
c
c                   END DO
c
c                 END DO
c              ENDDO

           ENDIF
        ENDDO
      ENDDO
#endif
C
C These two conditions pretty much take care of regular bonds. Now we
C need to consider isolated fragments and connect them.
C Identifying isolated fragments is a graph theory problem and I don't
C claim to be an expert on graph theory. There are efficient algorithms 
C and specialized programs to do this. The following is very simple
C and hopefully will work 90% of the time without heavy penalty. 
C 
C First, count the total number of connections for each atom and
C store the index of the other atoms that are connected to the
C given atom (in NCONPRCNTR, ICONNECTO arrays respectively). 
C Reorder these two arrays in descending order with respect to the
C number of connections.  
C Take a pair of atoms (starting with maximum connectivities) 
C and check whether they are connected. If they are not
C connected then search outwards (at predefined maxlevels
C representing a maximum possible chain length, for our
C cases I am thinking of 5). The searching is done as follows;
C Take the number of connections for the chosen pair and test whether at
C least one of them is connected; if so, then the choosen pair is not
C part of a disconnected fragment. If all of the atoms are not
C connected, then we look at the next level. When we consider all pairs
C and the FRAGMENT is still true then we have to mark that pair as a
C part of an isolated fragment. 
C 
      CALL IZERO(NCONPRCNTR, NRATMS)
      DO IBNDS = 1, NRATMS
         INDEX = 0
         DO JBNDS = 1, NRATMS
            IF (IBNDTO(JBNDS, IBNDS) .NE. JZERO .OR.
     &          IBNDTO(IBNDS, JBNDS) .NE. JZERO) THEN
                NCONPRCNTR(IBNDS) = NCONPRCNTR(IBNDS) + 1
                INDEX = INDEX + 1 
                ICONNECTO(IBNDS, INDEX) = JBNDS
                IF (NCONPRCNTR(IBNDS) .GT. MAXCNTVS) THEN
                    WRITE(6, 9) MAXCNTVS
                    CALL ERREX
                ENDIF
            ENDIF 
         ENDDO 
      ENDDO
  9   FORMAT(T4, 'The maximum no. of connectvities allowed for',
     &       'for a center ', I2, ' exceeded')
C
#ifdef _DEBUG_LVL
      Print*, "The number of connectivities per center"
      Print*,  (NCONPRCNTR(I), I=1, NRATMS)
      Print*, "The atom connection lables"
      Print*, ((ICONNECTO(I, J), J=1, NCONPRCNTR(I)),
     &         I=1, NRATMS)
#endif 
C

      do iatoms=2,nratms
         do jatoms=1,(iatoms-1)
            if (ibndto(jatoms,iatoms).eq.1) ibndto(iatoms,jatoms)=1

      call identfy_fragments(nratms,ibndto,connected_atoms,nfrags,
     &                       length_fragments,mark_fragments)
          enddo
      enddo
C
#ifdef _DEBUG_LVL0
      Print*, "The number of isolated fragments per center", nfrags
      Print*,  (length_fragments(I), I=1, NRATMS)       
      Print*, "Fragment atom lables"
      Print*, ((mark_fragments(I, J), J=1, nratms),
     &         I=1, NRATMS) 
      Print*, ((bndlenths(I, J), J=1, nratms),
     &         I=1, NRATMS) 
#endif
C
C Assign the primary interfragment connectivity
C
      IF (NFRAGS .GT. 1) THEN
         DO IATMS = 1, NFRAGS
            DO JATMS = NFRAGS, IATMS+1, -1
               ARCH_DIST = 10.0D0
               DO IFRAGS = 1, LENGTH_FRAGMENTS(IATMS)
                  DO JFRAGS = 1, LENGTH_FRAGMENTS(JATMS)
                     IFRAG_CENT  = MARK_FRAGMENTS(IATMS, IFRAGS)
                     JFRAG_CENT  = MARK_FRAGMENTS(JATMS, JFRAGS)
                     Print*, IFRAG_CENT, JFRAG_CENT
C
                      DIST_IJFRAG  = BNDLENTHS(IFRAG_CENT, 
     &                                         JFRAG_CENT)
                       Print*, "The distance", DIST_IJFRAG
                      IF (DIST_IJFRAG .LE. ARCH_DIST) THEN
                          ARCH_DIST = DIST_IJFRAG
                          I_ARCH    = IFRAG_CENT
                          J_ARCH    = JFRAG_CENT
                      ENDIF 
C                     BLNGT_IJFRG(JFRAGS) = DIST_IJFRAG
                   ENDDO
                ENDDO
#ifdef _DEBUG_LVL0
      Print*, "The smallest distance between fragments: ", ARCH_DIST
      Print*, "The atom indices of corrs. to smalest  : ", I_ARCH,
     &         J_ARCH 
#endif 
             IBNDTO(J_ARCH, I_ARCH) = IONE
             IBNDTO(I_ARCH, J_ARCH) = IONE 
             ENDDO
         ENDDO
      ENDIF
C
C Assign auxilary interfragment connectivities
C
      IF (NFRAGS .GT. 1) THEN
         DO IATMS = 1, NFRAGS
            DO JATMS = NFRAGS, IATMS+1, -1
               DO IFRAGS = 1, LENGTH_FRAGMENTS(IATMS)       
                  DO JFRAGS = 1, LENGTH_FRAGMENTS(JATMS)
                     IFRAG_CENT  = MARK_FRAGMENTS(IATMS, IFRAGS)
                     JFRAG_CENT  = MARK_FRAGMENTS(JATMS, JFRAGS)  
                     Print*, IFRAG_CENT, JFRAG_CENT
C
                     DIST_IJFRAG  = BNDLENTHS(IFRAG_CENT,
     &                                        JFRAG_CENT)
                     IF (DIST_IJFRAG .LT. 1.30D0*ARCH_DIST 
     &                   .OR. DIST_IJFRAG .LT. 2.0D0*
     &                   AMTOBOHR) THEN
                         IBNDTO(J_ARCH, I_ARCH) = IONE
                         IBNDTO(I_ARCH, J_ARCH) = IONE
                     ENDIF
                   ENDDO
                ENDDO
#ifdef _DEBUG_LVL0
      Print*, "The smallest distance between fragments: ", ARCH_DIST
      Print*, "The atom indices of corrs. to smalest  : ", I_ARCH,
     &         J_ARCH
#endif
             IBNDTO(J_ARCH, I_ARCH) = IONE
             IBNDTO(I_ARCH, J_ARCH) = IONE
             ENDDO
         ENDDO
      ENDIF
C
C Lets identify hydrogen bonds. H-bonds are assigned if H-X distance
C is within the sum of the covalent radii of H and X(=N, O, P, Cl, F,
C P, S) and 0.9 times the Van der Walls radii and Y...H-X angle
C is greater than 90 degrees.
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The connectivity assignments before H-bond block"
      Print*, ((ibndto(j, i), j=1, nratms), I=1, nratms)     
      Write(6,*)
#endif
      CALL DSCAL(NRATMS*NRATMS, 1.0D0/1.3D0, SMOFCOVRADI, 1)
      DO IBNDS = 1, NRATMS
         IF (IATOMICNMBER(IBNDS) .EQ. 1) THEN
            DO JBNDS = 1, (IBNDS - 1)
               NOPClFPS = ((IATOMICNMBER(JBNDS).EQ. 7)  .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 8)  .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 9)  .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 10) .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 15) .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 16) .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 17))

               IF (NOPClFPS.AND.IBNDTO(JBNDS,IBNDS).NE.JZERO) THEN
                  DO KBNDS = 1, NRATMS
                     IF (((IATOMICNMBER(KBNDS) .EQ. 7)  .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 8)  .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 9)  .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 10) .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 15) .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 16) .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 17)) .AND.
     &                   ((KBNDS .NE. IBNDS) .AND. (KBNDS .NE.
     &                    JBNDS))) THEN
                        IF (IBNDS .GT. KBNDS) THEN
                           DIST_HY     = BNDLENTHS(KBNDS, IBNDS)
                           SMOFCOVLRAD = SMOFCOVRADI(KBNDS,IBNDS)
                        ELSE
                           DIST_HY = BNDLENTHS(IBNDS, KBNDS)
                           SMOFCOVLRAD = SMOFCOVRADI(IBNDS,KBNDS)
                        END IF
                        CALL CMP_ANGLE(JBNDS, IBNDS, KBNDS,
     &                                 CARTCOORD, ANG_XHY, NRATMS)
                        KATOM = IATOMICNMBER(KBNDS)
                        IATOM = IATOMICNMBER(IBNDS)
                        SMOFVRADII = (VRADII(KATOM)+VRADII(IATOM))*0.9D0
#ifdef _DEBUG_LVL0
                        Print*, "Critical parameters for H-bonds"
                        Print*, DIST_HY, ANG_XHY , SMOFVRADII,
     &                          SMOFCOVLRAD
#endif
                        IF ((DIST_HY .GT. SMOFCOVLRAD).AND.
     &                      (DIST_HY .LT. SMOFVRADII).AND.
     &                      (ANG_XHY .GT. 90.0D0)) THEN
                           IF (IBNDS .GT. KBNDS) THEN
                              IBNDTO(KBNDS, IBNDS) = IONE
                           ELSE
                              IBNDTO(IBNDS, KBNDS) = IONE
                           END IF
                        END IF
                     END IF
                  END DO
               END IF
            END DO
         END IF
      END DO
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The connectivity assignments after H-bond block"
      Print*, ((ibndto(j, i), j=1, nratms),I=1, nratms)
      Write(6,*)
#endif
C 
C Symmetrize the connectivity array.
C
cSSS      DO 50 IATMS = 1, NRATMS
cSSS         DO 60 IBNDS = 1, NCONPRCNTR(IATMS)
cSSS            KBNDS = IBNDTO(IBNDS, IATMS)
cSSS            FLAG = .FALSE.
cSSS            DO 70 JBNDS = 1, NCONPRCNTR(KBNDS)
cSSS
cSSS                  IF (IBNDTO(JBNDS, KBNDS) .EQ. IATMS) THEN
cSSS                  FLAG = .TRUE.
cSSS                  CALL WRITE(LUOUT, 2000)
cSSS                  CALL ERREX
cSSS               ENDIF
cSSS
cSSS 70         CONTINUE
cSSS 60      CONTINUE
cSSS     
cSSS         IF (.NOT. FLAG) THEN
cSSS            NCONPRCNTR(KBNDS) = NCONPRCNTR(KBNDS) + 1
cSSS            IBNDTO(NCONPRCNTR(KBNDS), KBNDS) = IATMS
cSSS         ENDIF
cSSS
cSSS 50   CONTINUE
cSSS
cSSS 1000 FORMAT("@-assign_connectvities - Exceeds the maximum number of 
cSSS     &        bonds allowed for a center")
cSSS
cSSS 2000 FORMAT("@-assign_connectvitis - Duplicate connectivity", " I=",I5,
cSSS              " IB=",I3," J=",I5," JB=",I3,'.')

      RETURN
      END
 
