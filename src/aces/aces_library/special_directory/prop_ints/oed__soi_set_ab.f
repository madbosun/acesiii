         SUBROUTINE  OED__SOI_SET_AB
     +
     +                    ( NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      NUCLEI,
     +                      XN,YN,ZN,
     +                      IXDERC,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      DERCX,DERCY,DERCZ,
     +                      EXP1,EXP2,
     +                      CC1,CC2,
     +                      SPHERIC,
     +
     +                                 NCGTOA,NCGTOB,
     +                                 NPGTOA,NPGTOB,
     +                                 SHELLA,SHELLB,SHELLP,
     +                                 MXSHELL, ! need to be careful
     +                                 XA,YA,ZA,XB,YB,ZB,
     +                                 NDERX,NDERY,NDERZ,
     +                                 DERAX,DERAY,DERAZ,
     +                                 DERBX,DERBY,DERBZ,
     +                                 DIFFA,DIFFB,DIFFC,
     +                                 DIFFX,DIFFY,DIFFZ,
     +                                 IXAEQB,IXCEQA,IXCEQB,
     +                                 ATOMIC,EQUALAB,
     +                                 ABX,ABY,ABZ,RNABSQ,
     +                                 SPNORM,
     +                                 NXYZA,NXYZB,NXYZT,
     +                                 NRYA,NRYB,
     +                                 INDEXA,INDEXB,
     +               SWAP12,SWAPRS,   ! need to be careful
     +                                 LEXPA,LEXPB,
     +                                 LCCA,LCCB,
     +                                 LCCSEGA,LCCSEGB,
     +                                 EMPTY )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_SET_DERV_AB
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the logistics on how to evaluate
C                the (1|2) derivative nuclear attraction integral batch
C                in the most efficient way. It performs the label map:
C
C                               (1|2) --> (A|B)
C
C                and returns some additional data of crucial importance
C                for evaluation of the (A|B) derivative nuclear
C                attraction integrals.
C
C                The freedom we have in making the internal association
C                1,2 -> A,B follows from the 2-fold permutational
C                symmetry of the nuclear attraction integrals in (1|2):
C
C                                 (1|2) = (2|1)
C
C                Note, that this permutational symmetry only holds as
C                long as we permute the differential operators with it.
C
C                Application of the 1 <-> 2 switch is pretty irrelevant
C                here, due to the fact the routine which generates the
C                initial 1D integrals has the possibility to either
C                apply the HRR steps on the A or B shell. Still we
C                choose here to have A>=B for convenience and possibly
C                slightly better performance.
C
C
C                  Input (x = 1 and 2):
C
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                    NUCLEI       =  the original # of nuclear
C                                    attraction centers
C                    XN,YN,ZN     =  the original x,y,z-coordinates
C                                    of the nuclear attraction centers
C                    IXDERC       =  the index of which of the nuclear
C                                    attraction centers is to be
C                                    differentiated. If that index
C                                    corresponds to one of the centers
C                                    1 and/or 2, it will already be
C                                    differentiated along with these
C                                    centers, hence values transmitted
C                                    for DERCX,DERCY,DERCZ are
C                                    irrelevant in that case.
C                    DERyp        =  the order of differentiation on
C                                    centers y = 1 and 2 with respect
C                                    to the p = x,y,z coordinates
C                    DERCp        =  the order of differentiation for
C                                    the IXDERC-th nuclear attraction
C                                    center with respect to the
C                                    p = x,y,z coordinates
C                    EXPx         =  primitive exponents for csh x
C                    CCx          =  contraction coeffs for csh x
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C                  Output (x = A and B):
C
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    SHELLP       =  the shell sum A+B
C                    MXSHELL      =  the largest (maximum) shell type
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = A and B
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERyp        =  the order of differentiation on
C                                    centers y = A and B with respect
C                                    to the p = x,y,z coordinates
C                    DIFFy        =  is true, if differentiation will
C                                    be performed on centers y = A and B
C                                    and at least on one of the nuclear
C                                    attraction centers y = C different
C                                    from A and B and if differentiation
C                                    will be performed at least once
C                                    one on each of the y = x,y,z
C                                    coordinates
C                    IXAEQB       =  is an index = 0 or 1, depending
C                                    on if center A is different or
C                                    equal from center B, respectively
C                    IXCEQy       =  contains index of which of all
C                                    nuclear attraction centers is
C                                    equal to centers y = A,B. If no
C                                    such C is equal the value is 0
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated,
C                                    observing only centers A and B and
C                                    ignoring all nuclear attraction
C                                    centers
C                    EQUALAB      =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pair AB
C                    ABX,ABY,ABZ  =  the x,y,z-coordinate differences
C                                    between centers A and B
C                    RNABSQ       =  square of the magnitude of the
C                                    distance between centers A and B
C                    SPNORM       =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                    NXYZx        =  # of cartesian monomials for csh x
C                    NXYZT        =  total # cartesian monomial
C                                    combinations, that is NXYZA * NXYZB
C                    NRYx         =  # of spherical functions for csh x
C                    INDEXx       =  index A,B -> 1,2 map
C                    SWAP12       =  is .true., if a swap 1 <-> 2 has
C                                    been performed
C                    SWAPRS       =  is set .true. if the contraction
C                                    order of the primitives pair AB
C                                    will be performed in reverse order
C                                    BA for efficiency reasons
C                    LEXPx        =  pointers to locate appropriate
C                                    section of the exponent array
C                                    corresponding to csh x
C                    LCCx         =  pointers to locate appropriate
C                                    section of the contraction coeff
C                                    array corresponding to csh x
C                    LCCSEGx      =  pointers to locate appropriate
C                                    section of the lowest and highest
C                                    primitive index array defining
C                                    segmented contraction boundaries
C                                    for csh x
C                    EMPTY        =  logical flag, indicating if an
C                                    empty batch of integrals is
C                                    expected.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ALERT
         LOGICAL     ATOMIC,ATOMAC,ATOMBC
         LOGICAL     DIFFA,DIFFB,DIFFC
         LOGICAL     DIFFX,DIFFY,DIFFZ
         LOGICAL     EMPTY
         LOGICAL     EQUALAB
         LOGICAL     NODIFF
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAPRS

         INTEGER     DER1X,DER1Y,DER1Z
         INTEGER     DER2X,DER2Y,DER2Z
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DERCX,DERCY,DERCZ
         INTEGER     I,J,N
         INTEGER     INDEXA,INDEXB
         INTEGER     IXAEQB,IXCEQA,IXCEQB
         INTEGER     IXDERC
         INTEGER     LCCSEGA,LCCSEGB
         INTEGER     LCCA,LCCB
         INTEGER     LEXPA,LEXPB
         INTEGER     MXSHELL
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NCGTOA,NCGTOB
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NRYA,NRYB
         INTEGER     NUCLEI
         INTEGER     NXYZA,NXYZB,NXYZT
         INTEGER     SHELL1,SHELL2
         INTEGER     SHELLA,SHELLB
         INTEGER     SHELLP

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  RNABSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  XN (1:NUCLEI)
         DOUBLE PRECISION  YN (1:NUCLEI)
         DOUBLE PRECISION  ZN (1:NUCLEI)

         DOUBLE PRECISION  EXP1 (1:NPGTO1)
         DOUBLE PRECISION  EXP2 (1:NPGTO2)

         DOUBLE PRECISION  CC1 (1:NPGTO1,1:NCGTO1)
         DOUBLE PRECISION  CC2 (1:NPGTO2,1:NCGTO2)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...generate all 1,2 data.
C
C
         EMPTY  = .FALSE.
         ATOMIC = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         SHELLP = SHELL1 + SHELL2

c         write(*,*) 'X1 and X2', X1, X2, (X1.EQ.X2)
c         write(*,*) 'Y1 and Y2', Y1, Y2, (Y1.EQ.Y2)
c         write(*,*) 'Z1 and Z2', Z1, Z2, (Z1.EQ.Z2)

c      write(*,*) 'Shellp= she1 + shel2',SHELLP,shell1, shell2
         MXSHELL = MAX0 (SHELL1,SHELL2)
c      write(*,*) 'MXSHELL=',MXSHELL

cPV I do not think we will need MXSHELL info 

C
C
C             ...determine csh equality between centers 1 and 2
C                in increasing order of complexity:
C
C                 centers -> shells -> exponents -> ctr coefficients
C
C
         EQUALAB = ATOMIC

         IF (EQUALAB) THEN
             EQUALAB =     (SHELL1 .EQ. SHELL2)
     +               .AND. (NPGTO1 .EQ. NPGTO2)
     +               .AND. (NCGTO1 .EQ. NCGTO2)
             IF (EQUALAB) THEN
               DO I = 1,NPGTO1
                  EQUALAB = EQUALAB .AND. (EXP1(I).EQ.EXP2(I))
               END DO
               IF (EQUALAB) THEN
                 DO J = 1,NCGTO1
                    IF (EQUALAB) THEN
                      DO I = 1,NPGTO1
                         EQUALAB = EQUALAB .AND. (CC1(I,J).EQ.CC2(I,J))
                      END DO
                    END IF
                 END DO
               END IF
             END IF
         END IF
CSSS       write(*,*) 'force equlab to be false'
        EQUALAB=.false.
C
C
C             ...decide on the 1 <-> 2 swapping.
C
C
         SWAP12 = SHELL1 .LT. SHELL2
         swap12=.false.
cPV SWAP12 is redundant I can not use SWAP12 logic


C
C
C             ...according to the previously gathered info, set the
C                new A,B shells, # of primitives + contraction coeffs
C                as well as pointers to the alpha exponents and
C                contraction coefficients.
C
C
         IF (.NOT.SWAP12) THEN
             XA = X1
             YA = Y1
             ZA = Z1
             XB = X2
             YB = Y2
             ZB = Z2
             SHELLA = SHELL1
             SHELLB = SHELL2
             NPGTOA = NPGTO1
             NPGTOB = NPGTO2
             NCGTOA = NCGTO1
             NCGTOB = NCGTO2
             DERAX = DER1X
             DERAY = DER1Y
             DERAZ = DER1Z
             DERBX = DER2X
             DERBY = DER2Y
             DERBZ = DER2Z
             INDEXA = 1
             INDEXB = 2
             LEXPA = 1
             LEXPB = LEXPA + NPGTO1
             LCCA = 1
             LCCB = LCCA + NPGTO1 * NCGTO1
             LCCSEGA = 1
             LCCSEGB = LCCSEGA + NCGTO1
         ELSE
cPV I need to ignore this block
             XA = X2
             YA = Y2
             ZA = Z2
             XB = X1
             YB = Y1
             ZB = Z1
             SHELLA = SHELL2
             SHELLB = SHELL1
             NPGTOA = NPGTO2
             NPGTOB = NPGTO1
             NCGTOA = NCGTO2
             NCGTOB = NCGTO1
             DERAX = DER2X
             DERAY = DER2Y
             DERAZ = DER2Z
             DERBX = DER1X
             DERBY = DER1Y
             DERBZ = DER1Z
             INDEXA = 2
             INDEXB = 1
             LEXPB = 1
             LEXPA = LEXPB + NPGTO1
             LCCB = 1
             LCCA = LCCB + NPGTO1 * NCGTO1
             LCCSEGB = 1
             LCCSEGA = LCCSEGB + NCGTO1
cPV end the ignore block
         END IF
C
C
C             ...the new A and B shells are set. Calculate the
C                following info:
C
C                1) control variable to determine contraction order
C                2) cartesian monomial dimensions
C                3) spherical dimensions (= cartesian, if no spherical)
C                4) the overall norm factor due to s- or p-type shells
C
C

cPV SWAPRS logic variable is redunant as well


         SWAPRS = NPGTOA .GT. NPGTOB
        
CSSS         write(*,*) 'SWAPRS=',SWAPRS
cPV end

         NXYZA  = (SHELLA+1)*(SHELLA+2)/2
         NXYZB  = (SHELLB+1)*(SHELLB+2)/2
         NXYZT  = NXYZA * NXYZB

         NRYA = SHELLA + SHELLA + 1
         NRYB = SHELLB + SHELLB + 1

         IF (.NOT.SPHERIC) THEN
             NRYA = NXYZA
             NRYB = NXYZB
         END IF

         SPNORM = ONE
         IF (SHELLA.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELLB.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
C
C
C             ...calculate the coordinate differences between centers
C                A and B and calculate the square of the magnitude
C                of the distance between A and B. Set also the A,B
C                equality index and determine if differentiation is
C                is to be performed on centers A and/or B.
C
C                Check here also, if the specified derivative operator
C                orders for A and B match the center equality patterns
C                for each x,y,z-component. Stop the program if any
C                inconsistency is found.
C
C

CSSS         Write(*,*) 'ATOMIC=',ATOMIC
         IF (.NOT.ATOMIC) THEN

             ABX = XA - XB
             ABY = YA - YB
             ABZ = ZA - ZB
             RNABSQ = ABX * ABX + ABY * ABY + ABZ * ABZ
             IXAEQB = 0

             DIFFA = (DERAX + DERAY + DERAZ) .NE. 0
             DIFFB = (DERBX + DERBY + DERBZ) .NE. 0

             NDERX = DERAX + DERBX
             NDERY = DERAY + DERBY
             NDERZ = DERAZ + DERBZ

         ELSE

             ABX = ZERO
             ABY = ZERO
             ABZ = ZERO
             RNABSQ = ZERO
             IXAEQB = 1

CSSS             ALERT =     DERAX.NE.DERBX
CSSS     +              .OR. DERAY.NE.DERBY
CSSS     +              .OR. DERAZ.NE.DERBZ

c             Write(*,'(a,6(1x,I2))') 'Prakash DERAX,DERBX,
c     +                    DERAYR, DERBY,
c     +                    DERAZ, DERBZ',
c     +                    DERAX,DERBX,
c     +                    DERAY,DERBY,
c     +                    DERAZ,DERBZ
CSSS             IF (ALERT) THEN
CSSS                 WRITE (*,*) ' Center A = B / derv order wrong! '
CSSS                 WRITE (*,*) ' oed__nai_set_derv_ab '
cPV                 STOP
CSSS             END IF

             DIFFA = (DERAX + DERAY + DERAZ) .NE. 0
cPV             DIFFB = DIFFA
             DIFFB=(DERBX + DERBY + DERBZ) .NE. 0

cPV             NDERX = DERAX
cPV             NDERY = DERAY
cPV             NDERZ = DERAZ

             NDERX = DERAX + DERBX
             NDERY = DERAY + DERBY
             NDERZ = DERAZ + DERBZ
        
         END IF

         DIFFX = (DERAX + DERBX) .NE. 0
         DIFFY = (DERAY + DERBY) .NE. 0
         DIFFZ = (DERAZ + DERBZ) .NE. 0
c         Write(*,*)'DIFFX=',DIFFx,diffy,diffz

C
C
C             ...complete next the differentiation info due to the
C                possible differentiation of a nuclear attraction
C                center C:
C
C                   i) x,y,z-component and C center differentiation
C                      indicator
C
C                  ii) the overall differentiation order, observing
C                      possible center equalities. Nuclear attraction
C                      centers C may be equal to centers A and/or B
C                      but they are always mutually different among
C                      each other.
C
C                 iii) index of nuclear attraction center C equal to
C                      A and/or B. If no center C is equal to either
C                      A and/or B then this index is zero. If more
C                      than one nuclear attraction center C is equal to
C                      A and/or B the program stops with a diagnostic.
C
C
         IXCEQA = 0
         IXCEQB = 0

         DO N = 1,NUCLEI

            XC = XN (N)
            YC = YN (N)
            ZC = ZN (N)

            ATOMAC = (XA.EQ.XC) .AND. (YA.EQ.YC) .AND. (ZA.EQ.ZC)
            ATOMBC = (XB.EQ.XC) .AND. (YB.EQ.YC) .AND. (ZB.EQ.ZC)

            IF (ATOMAC) THEN
                IF (IXCEQA.NE.0) THEN
                    WRITE (*,*) ' Nuclear attraction centers equal! '
                    WRITE (*,*) ' oed__nai_set_derv_ab '
                    STOP
                END IF
                IXCEQA = N
            END IF

            IF (ATOMBC) THEN
                IF (IXCEQB.NE.0) THEN
                    WRITE (*,*) ' Nuclear attraction centers equal! '
                    WRITE (*,*) ' oed__nai_set_derv_ab '
                    STOP
                END IF
                IXCEQB = N
            END IF

         END DO

         DIFFC = (DERCX+DERCY+DERCZ) .NE. 0
     +           .AND. (IXDERC.NE.IXCEQA)
     +           .AND. (IXDERC.NE.IXCEQB)
     +           .AND. (IXDERC.GT.0)
     +           .AND. (IXDERC.LE.NUCLEI)

         IF (DIFFC) THEN
             NDERX = NDERX + DERCX
             NDERY = NDERY + DERCY
             NDERZ = NDERZ + DERCZ
             DIFFX = DIFFX .OR. (DERCX.NE.0)
             DIFFY = DIFFY .OR. (DERCY.NE.0)
             DIFFZ = DIFFZ .OR. (DERCZ.NE.0)
         END IF

         NODIFF = (NDERX.EQ.0) .AND. (NDERY.EQ.0) .AND. (NDERZ.EQ.0)

         IF (NODIFF) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
c         write(*,*) 'Diffc=',diffC,diffx,diffy,diffz
C
C
C             ...ready!
C
C
         RETURN
         END
