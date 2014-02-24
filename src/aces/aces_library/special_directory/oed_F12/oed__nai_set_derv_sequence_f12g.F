         SUBROUTINE  OED__NAI_SET_DERV_SEQUENCE_F12G
     +
     +                    ( NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DERCX,DERCY,DERCZ,
     +                      DIFFA,DIFFB,DIFFC,
     +                      ATOMIC,
     +
     +                               CENSQX,CENSQY,CENSQZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_SET_DERV_SEQUENCE
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine determines the sequence of centers for
C                the x,y,z-coordinate derivations. The info is stored
C                in arrays CENSQp for each coordinate p=X,Y,Z.
C
C                The indices for the centers obey the following
C                indexing scheme:
C
C                               centers: A  B  C
C                          center index: 1  2  3
C
C                where C denotes the extra nuclear attraction center
C                different from A and B that is to be differentiated.
C
C                Take the x-coordinate for example and let CENSQX have
C                the following elements on exit from this routine:
C
C                            CENSQX (I) = 2,3,3,1,1
C
C                Then the sequence of single derivations is:
C
C                   dx/dAx  dx/dAx  dx/dCx  dx/dCx  dx/dBx  (A|B)
C
C                Note, that the indexing of the CENSQp arrays start
C                with 0. The reason is that the NDERp values might
C                be equal to zero.
C
C                The strategy to determine the center sequences is
C                to perform differentiation on the centers A and/or B
C                first, since they reduce the angular momentum by -1
C                for the 1D integrals involved. On the other hand,
C                differentiation on any of the nuclear attraction
C                centers not equal to A and B does not change the
C                1D angular momentum parts on centers A and B.
C
C
C                  Input:
C
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERxp        =  the order of differentiation on
C                                    centers x = A and B with respect
C                                    to the p = x,y,z coordinates
C                    DERCp        =  the order of differentiation for
C                                    the IXDERC-th nuclear attraction
C                                    center with respect to the
C                                    p = x,y,z coordinates
C                    DIFFy        =  is true, if differentiation will
C                                    be performed on centers y = A and B
C                                    and at least on one of the nuclear
C                                    attraction centers y = C different
C                                    from A and B
C                    ATOMIC       =  indicates, if centers A and B are
C                                    equal
C
C                  Output:
C
C                    CENSQp       =  center sequence array for the
C                                    p = x,y,z coordinates
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C                      >>>  COPYRIGHT NOTICE  <<<
C                             Program ERD
C
C    This routine is part of the integral program package ERD and as
C    such is copyrighted by N. Flocke, 2004; all rights reserved.
C    The full copyright notice is printed in the main routine of the
C    ERD package and all terms and conditions given therein are
C    incorporated here.
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC
         LOGICAL     DIFFA,DIFFB,DIFFC

         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DERCX,DERCY,DERCZ
         INTEGER     I
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NX,NY,NZ

         INTEGER     CENSQX (0:NDERX)
         INTEGER     CENSQY (0:NDERY)
         INTEGER     CENSQZ (0:NDERZ)
C
C
C------------------------------------------------------------------------
C
C
C             ...first the center sequences due to centers A.
C
C
         NX = 0
         NY = 0
         NZ = 0

         IF (DIFFA) THEN

             DO I = 1,DERAX
                CENSQX (NX+I) = 1
             END DO
             NX = NX + DERAX

             DO I = 1,DERAY
                CENSQY (NY+I) = 1
             END DO
             NY = NY + DERAY

             DO I = 1,DERAZ
                CENSQZ (NZ+I) = 1
             END DO
             NZ = NZ + DERAZ

         END IF
C
C
C             ...add the center sequences due to centers B (if needed).
C
C
         IF (DIFFB .AND. .NOT.ATOMIC) THEN

             DO I = 1,DERBX
                CENSQX (NX+I) = 2
             END DO
             NX = NX + DERBX

             DO I = 1,DERBY
                CENSQY (NY+I) = 2
             END DO
             NY = NY + DERBY

             DO I = 1,DERBZ
                CENSQZ (NZ+I) = 2
             END DO
             NZ = NZ + DERBZ

         END IF
C
C
C             ...center sequences due to a possible differentiation
C                on a nuclear attraction center different from centers
C                A and B.
C
C
         IF (DIFFC) THEN

             DO I = 1,DERCX
                CENSQX (NX+I) = 3
             END DO
             NX = NX + DERCX

             DO I = 1,DERCY
                CENSQY (NY+I) = 3
             END DO
             NY = NY + DERCY

             DO I = 1,DERCZ
                CENSQZ (NZ+I) = 3
             END DO
             NZ = NZ + DERCZ

         END IF
C
C
C             ...safety check.
C
C
         IF (NX.NE.NDERX .OR. NY.NE.NDERY . OR. NZ.NE.NDERZ) THEN
             WRITE (*,*) ' Total derivative order mismatch! '
             WRITE (*,*) ' NX,NY,NZ,NDERX,NDERY,NDERZ = ',
     +                     NX,NY,NZ,NDERX,NDERY,NDERZ
             WRITE (*,*) ' oed__nai_set_derv_sequence '
             STOP
         END IF

C
C             ...ready!
C
C
         RETURN
         END
