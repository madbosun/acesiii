         SUBROUTINE  OED__NAI_DERV_PCGTO_BLOCK_F12G
     +
     +                    ( NBATCH,
     +                      NINT1DX,NINT1DY,NINT1DZ,
     +                      ATOMIC,
     +                      MIJ,NCEN,MIJCEN,
     +                      NIJ,NIJBEG,NIJEND,
     +                      NGQP,NMOM,NGQSCR,MGIJCEN,
     +                      NPGTOA,NPGTOB,
     +                      NXYZA,NXYZB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      NUCLEI,
     +                      XN,YN,ZN,NCHARGE,
     +                      NUCCEN,
     +                      NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      DIFFA,DIFFB,DIFFC,
     +                      IXAEQB,IXCEQA,IXCEQB,
     +                      ALPHAA,ALPHAB,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      CENSQX,CENSQY,CENSQZ,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      PX,PY,PZ,PAX,PAY,PAZ,PINVHF,SCALE,
     +                      RTS,WTS,GQSCR,TVAL,UVAL,
     +                      R1X,R1Y,R1Z,R2,
     +                      EXP2A,EXP2B,EXP2AB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      INT1DX_2,INT1DY_2,INT1DZ_2,STGEXP,
     +                      RYSGEN,CORETOP,RYSDATA,RYSDATA1,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__RYS_ROOTS_WEIGHTS
C                OED__NAI_1D_COEFFICIENTS
C                OED__NAI_1D_AB_INTEGRALS
C                OED__NAI_1D_DERV_INTEGRALS
C                OED__NAI_DERV_INT1D_TO_00
C                OED__NAI_DERV_INT1D_TO_A0
C                OED__NAI_DERV_INT1D_TO_AB
C  DESCRIPTION : This operation calculates a batch of derivated unnormed
C                nuclear attraction integrals between cartesian
C                gaussians
C
C                                   [A|B]
C                                        ij
C
C                for a block of ij exponent pairs. The total number
C                of nuclear attraction integrals generated here is thus
C                given by the total number of cartesian monomials
C                NXYZA * NXYZB times the total number of exponent pairs
C                MIJ in the present block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,nxyza*nxyzb)
C
C
C                  Input:
C
C                    NBATCH       =  size of the array that will hold
C                                    the final primitive cartesian
C                                    derivative integral batch as well
C                                    as intermediate differentiated
C                                    1D integrals
C                    NINT1Dx      =  space needed for each of the 1D
C                                    x = X,Y,Z integral arrays (they
C                                    might be different due to different
C                                    orders of differentiation for
C                                    each cartesian component)
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
C                    NCEN         =  # of nuclear attraction centers
C                    MIJCEN       =  # of ij primitive index pairs
C                                    times # of nuclear attraction
C                                    centers
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    NGQSCR       =  size of gaussian quadrature
C                                    scratch space needed to calculate
C                                    all the quadrature roots
C                    MGIJCEN      =  # of roots times # of ij primitive
C                                    index pairs times # of nuclear
C                                    attraction centers
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    NXYZx        =  # of cartesian monomials for
C                                    each contraction shell x = A,B
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A,B and P=A+B
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B
C                    ABm          =  the m=x,y,z-coordinate differences
C                                    between centers A and B
C                    NUCLEI       =  # of nuclear attraction centers
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    NCHARGE      =  the nuclear charges for all
C                                    nuclear attraction centers
C                    NUCCEN       =  contains those index labels of
C                                    the nuclear attraction centers
C                                    that survived the screening process
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
C                    IXCEQy       =  contains the index of which
C                                    screened nuclear attraction center
C                                    C is equal to centers y = A,B.
C                                    If no C is equal, the value is 0
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    CENSQp       =  derivative center sequence array
C                                    for the p = x,y,z coordinates
C                    PRIMx        =  i,j labels of primitives for the
C                                    respective contraction shells
C                                    x = A,B
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B
C                    RHOAB        =  the complete set of NIJ exponential
C                                    prefactors between contraction
C                                    shells A and B
C                    Px           =  will hold current MIJ coordinates
C                                    x=X,Y,Z for the gaussian product
C                                    centers P=A+B
C                    PAx          =  will hold current MIJ coordinate
C                                    x=X,Y,Z differences P-A between
C                                    centers P and A
C                    PINVHF       =  will hold current MIJ values of
C                                    1/(2*P), where P are the exponent
C                                    sums for contraction shells A
C                                    and B
C                    SCALE        =  will hold current distinct MIJCEN
C                                    (expanded to MGIJCEN) values of
C                                    scaling factors
C                    RTS          =  will hold all current MGIJCEN
C                                    quadrature roots
C                    WTS          =  will hold all current MGIJCEN
C                                    quadrature weights
C                    GQSCR        =  will be used as scratch space
C                                    for determining the quadrature
C                                    roots and weights
C                    TVAL         =  will hold current MIJCEN values
C                                    of T-exponents defining the Rys
C                                    weight functions
C                    R1x          =  will hold the current MGIJCEN
C                                    VRR R1-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center P
C                    R2           =  will hold the current MGIJCEN
C                                    coordinate independent VRR
C                                    R2-coefficients
C                    EXP2x        =  will hold current double primitive
C                                    exponent values in MIJ order
C                                    (expanded to MGIJCEN) for each
C                                    contracted shell x = A,B
C                    EXP2AB       =  will hold current double primitive
C                                    exponent sum values in MIJ order
C                                    between both contracted shells
C                                    A and B
C                    INT1Dx       =  will hold all current derivated
C                                    1D integrals for each cartesian
C                                    component (x = X,Y,Z) during all
C                                    stages of differentiation
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    derivative cartesian nuclear
C                                    attraction [A|B] integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C                      >>>  COPYRIGHT NOTICE  <<<
C                             Program OED
C
C    This routine is part of the integral program package OED and as
C    such is copyrighted by N. Flocke, 2004; all rights reserved.
C    The full copyright notice is printed in the main routines of the
C    OED package and all terms and conditions given therein are
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
         LOGICAL     CEQA,CEQB
         LOGICAL     DIFFA,DIFFB,DIFFC
         LOGICAL     DIFFX,DIFFY,DIFFZ

         INTEGER     CENTER
         INTEGER     DERA,DERB
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DEROPA,DEROPB,DEROPC
         INTEGER     G000,G010,G020,G030,G040,G050,G060
         INTEGER     I,J,L,M,N
         INTEGER     IJ
         INTEGER     IXC,IXAEQB,IXCEQA,IXCEQB
         INTEGER     MGRID,NGRID
         INTEGER     MIJ,MGQPIJ,MIJCEN,MGIJCEN
         INTEGER     NBATCH
         INTEGER     NC,NCEN
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NG,NGQP,NGQPCEN,NMOM,NGQSCR
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NINT,NINT1DX,NINT1DY,NINT1DZ
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NUCLEI,RYSGEN,CORETOP
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB,SHELLP

         INTEGER     CENSQX (0:NDERX)
         INTEGER     CENSQY (0:NDERY)
         INTEGER     CENSQZ (0:NDERZ)

         INTEGER     NUCCEN (1:NCEN)

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  CX,CY,CZ
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  P,PINV,PVAL
         DOUBLE PRECISION  PCX,PCY,PCZ
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL
         DOUBLE PRECISION  RNPCSQ
         DOUBLE PRECISION  SCALEM,SCALEN
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  ZERO,HALF,ONE

         DOUBLE PRECISION  XN      (1:NUCLEI)
         DOUBLE PRECISION  YN      (1:NUCLEI)
         DOUBLE PRECISION  ZN      (1:NUCLEI)
         DOUBLE PRECISION  NCHARGE (1:NUCLEI)

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  EXP2AB  (1:MIJ)
         DOUBLE PRECISION  PX      (1:MIJ)
         DOUBLE PRECISION  PY      (1:MIJ)
         DOUBLE PRECISION  PZ      (1:MIJ)
         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)

         DOUBLE PRECISION  GQSCR   (1:NGQSCR)
         DOUBLE PRECISION  TVAL    (1:MIJCEN)
         DOUBLE PRECISION  UVAL    (1:MIJCEN)

         DOUBLE PRECISION  R1X     (1:MGIJCEN)
         DOUBLE PRECISION  R1Y     (1:MGIJCEN)
         DOUBLE PRECISION  R1Z     (1:MGIJCEN)
         DOUBLE PRECISION  R2      (1:MGIJCEN)
         DOUBLE PRECISION  EXP2A   (1:MGIJCEN)
         DOUBLE PRECISION  EXP2B   (1:MGIJCEN)
         DOUBLE PRECISION  RTS     (1:MGIJCEN)
         DOUBLE PRECISION  SCALE   (1:MGIJCEN)
         DOUBLE PRECISION  WTS     (1:MGIJCEN)

         DOUBLE PRECISION  INT1DX  (1:NINT1DX)
         DOUBLE PRECISION  INT1DY  (1:NINT1DY)
         DOUBLE PRECISION  INT1DZ  (1:NINT1DZ)

         DOUBLE PRECISION  INT1DX_2  (1:NINT1DX)
         DOUBLE PRECISION  INT1DY_2  (1:NINT1DY)
         DOUBLE PRECISION  INT1DZ_2  (1:NINT1DZ)
         DOUBLE PRECISION  RYSDATA(*),RYSDATA1(*)


         DOUBLE PRECISION  FTABLE  (0:MGRID,0:NGRID),STGEXP

         PARAMETER  (ZERO = 0.D0)
         PARAMETER  (HALF = 0.5D0)
         PARAMETER  (ONE  = 1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the quantities needed to establish the
C                1D nuclear attraction integrals.
C
C
         IF (ATOMIC) THEN
             L = 0
             M = 0
             N = 0
             DO IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                P = ALPHAA (I) + ALPHAB (J)
                PX (M) = XA
                PY (M) = YA
                PZ (M) = ZA
                PINV = ONE / P
                PINVHF (M) = HALF * PINV
                SCALEM = PINV * NORMA (I) * NORMB (J)
                DO NC = 1,NCEN
                   N = N + 1
                   IXC = NUCCEN (NC)
                   PCX = XA - XN (IXC)
                   PCY = YA - YN (IXC)
                   PCZ = ZA - ZN (IXC)
                   RNPCSQ = PCX * PCX + PCY * PCY + PCZ * PCZ
                   TVAL (N) = P * RNPCSQ
                   UVAL (N) = (STGEXP*STGEXP)/(4.d0*P)
                   SCALEN = SCALEM * NCHARGE (IXC)
                   DO NG = 1,NGQP
                      L = L + 1
                      SCALE (L) = SCALEN
                   END DO
                END DO
             END DO
         ELSE
             L = 0
             M = 0
             N = 0
             DO IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                EXPA = ALPHAA (I)
                EXPB = ALPHAB (J)
                P = EXPA + EXPB
                PINV = ONE / P
                PVAL = - EXPB * PINV
                PAX (M) = PVAL * ABX
                PAY (M) = PVAL * ABY
                PAZ (M) = PVAL * ABZ
                PXVAL = PAX (M) + XA
                PYVAL = PAY (M) + YA
                PZVAL = PAZ (M) + ZA
                PX (M) = PXVAL
                PY (M) = PYVAL
                PZ (M) = PZVAL
                PINVHF (M) = HALF * PINV
                SCALEM = PINV * NORMA (I) * NORMB (J) * RHOAB (IJ)
                DO NC = 1,NCEN
                   N = N + 1
                   IXC = NUCCEN (NC)
                   PCX = PXVAL - XN (IXC)
                   PCY = PYVAL - YN (IXC)
                   PCZ = PZVAL - ZN (IXC)
                   RNPCSQ = PCX * PCX + PCY * PCY + PCZ * PCZ
                   TVAL (N) = P * RNPCSQ
                   UVAL (N) = (STGEXP*STGEXP)/(4.d0*P)
                   SCALEN = SCALEM * NCHARGE (IXC)
                   DO NG = 1,NGQP
                      L = L + 1
                      SCALE (L) = SCALEN
                   END DO
                END DO
             END DO
         END IF
C
C
C             ...determine memory allocation offsets for the scratch
C                arrays used to calculate the quadrature roots +
C                weights:
C
C                   G000 = offset for A coefficients (Jacobi/Laguerre)
C                   G010 = offset for B coefficients (Jacobi/Laguerre)
C                   G020 = offset for moments (Jacobi/Laguerre)
C                   G030 = offset for diagonals of symmetric termat
C                   G040 = offset for offdiagonals of symmetric termat
C                   G050 = offset for first row intermediates during
C                          evaluation of symmetric termat
C                   G060 = offset for second row intermediates during
C                          evaluation of symmetric termat
C
C
         G000 = 1
         G010 = G000 + NMOM
         G020 = G010 + NMOM - 1
         G030 = G020 + NMOM
         G040 = G030 + NGQP
         G050 = G040 + NGQP
         G060 = G050 + NMOM
c*****************************************************************************
c   Here we do the most difficult part of the computation of F12             |
c   potential integrals - calculation of the roots and weight of the         |
c  corresponding Rys polynomials. The calculation of the necessary roots and |
c   weight is done by 2D ineterpolation (implemented by Toru Shiozaki)       |
c*****************************************************************************
C
C
C             ...calculate all roots and weights for all ij pairs and
C                nuclear centers. Array R2 is passed as a scratch array.
C
C

                  IF(RYSGEN.eq.1) then      
             

         CALL    OED__RYS_ROOTS_WEIGHTS_F12
     +
     +                ( MIJCEN,MGIJCEN,
     +                  NGQP,NMOM,
     +                  TVAL,UVAL,R2,
     +                  FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                  GQSCR(G000),GQSCR(G010),
     +                  GQSCR(G020),
     +                  GQSCR(G030),GQSCR(G040),
     +                  GQSCR(G050),GQSCR(G060),
     +
     +                           RTS,
     +                           WTS )
     +
     +

          DO I=1,MGIJCEN
          RYSDATA(I) = RTS(I)
          RYSDATA1(I) = WTS(I)
          ENDDO

                   ELSE

          DO I=1,MGIJCEN
          RTS(I) = RYSDATA(I)
          WTS(I) = RYSDATA1(I)
          ENDDO

                  ENDIF


c--------------------------------------------------------------------------------



C--------------------------------------------------------------------------------
C             ...calculate the derivative 1D integral coefficents               |
C                (if needed). These are the double exponent values              |
C                for differentiation on centers A and/or B and the              |
C                double exponent sum values between centers A and B.            |
C--------------------------------------------------------------------------------
         CEQA = IXCEQA .GT. 0
         CEQB = IXCEQB .GT. 0

         MGQPIJ = NGQP * MIJ
         NGQPCEN = NGQP * NCEN

         IF (DIFFA) THEN
             M = 0
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                EXPA = ALPHAA (I) + ALPHAA (I)
                DO N = 1,NGQPCEN
                   EXP2A (M+N) = EXPA
                END DO
                M = M + NGQPCEN
             END DO
         END IF

         IF (DIFFB) THEN
             M = 0
             DO IJ = 1,MIJ
                J = PRIMB (IJ)
                EXPB = ALPHAB (J) + ALPHAB (J)
                DO N = 1,NGQPCEN
                   EXP2B (M+N) = EXPB
                END DO
                M = M + NGQPCEN
             END DO
         END IF

         IF (DIFFC .OR. CEQA .OR. CEQB) THEN
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                J = PRIMB (IJ)
                P = ALPHAA (I) + ALPHAB (J)
                EXP2AB (IJ) = P + P
             END DO
         END IF
C----------------------------------------------------------------|
C                                                                |
C             ...generate all 1D VRR coefficients.               |
C                                                                |
C----------------------------------------------------------------|
         CALL    OED__NAI_1D_COEFFICIENTS
     +
     +                ( NGQP,NCEN,
     +                  MIJ,MGIJCEN,
     +                  ATOMIC,
     +                  NUCLEI,
     +                  XN,YN,ZN,
     +                  NUCCEN,
     +                  PX,PY,PZ,
     +                  PAX,PAY,PAZ,
     +                  PINVHF,
     +                  RTS,
     +
     +                            R1X,R1Y,R1Z,
     +                            R2 )
     +
     +
C-------------------------------------------------------------------------
C                                                                        |
C             ...start assembling the 1D AB nuclear attraction           |
C                integrals and their differentiation. Assemble first     |
C                the 1DX AB integrals and perform differentiation        |
C                sequence, if necessary.                                 |
C                                                                        |
C-------------------------------------------------------------------------
         IF (DIFFX) THEN

            DERA = DERAX
            DERB = DERBX



             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1X,R2,
     +                      ABX,
     +                      .TRUE.,
     +                      BATCH,
     +
     +                                INT1DX )

             CENSQX (1)=3
             CENSQX (2)=1

             DO N = 1,NDERX


                NINT = MGIJCEN * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DX (M)
                END DO

                CENTER = CENSQX (N)
c                 write(6,*) 'CENTER =',CENTER 
c                  IF(M.eq.1) CENTER=3
c                  IF(M.eq.2) CENTER=1 

                IF (CENTER.EQ.1) THEN

                    DEROPA = 1
                    DEROPB = IXAEQB
                    DEROPC = 0
                    IF (CEQA) THEN
                        DEROPC = 1
                        IXC = IXCEQA
                        CX = XN (NUCCEN (IXCEQA))
                    END IF

                ELSE IF (CENTER.EQ.2) THEN

                    DEROPA = IXAEQB
                    DEROPB = 1
                    DEROPC = 0
                    IF (CEQB) THEN
                        DEROPC = 1
                        IXC = IXCEQB
                        CX = XN (NUCCEN (IXCEQB))
                    END IF

                ELSE

                    DEROPA = 0
                    DEROPB = 0
                    DEROPC = 1
                    IXC = 1
                    CX = XN (NUCCEN (1))

                END IF

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB


                IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                    WRITE (*,*) ' Problems forming x-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__nai_derv_pcgto_block '
                    STOP
                END IF

c                write(6,*) DERA,DERB,DEROPC
c                write(6,*) '<========='

                CALL    OED__NAI_1D_DERV_INTEGRALS
     +
     +                       ( MIJ,MGQPIJ,MGIJCEN,
     +                         NGQP,NGQPCEN,
     +                         SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,DEROPC,
     +                         EXP2A,EXP2B,EXP2AB,
     +                         IXC,PX,CX,
     +                         RTS,
     +                         R1X,
     +                         BATCH,
     +
     +                                  INT1DX )
     +
     +
             END DO
c            write(6,*) '----------------------------------------'
c======================================================================



c------------------------------------------------------------------------------

         ELSE
             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1X,R2,
     +                      ABX,
     +                      .TRUE.,
     +                      BATCH,
     +
     +                                INT1DX )
     +
     +
         END IF
C--------------------------------------------------------------
C                                                             |
C             ...assemble the 1DY AB integrals and perform    |
C                differentiation sequence, if necessary.      |
C                                                             |
C--------------------------------------------------------------
         IF (DIFFY) THEN

             DERA = DERAY
             DERB = DERBY

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Y,R2,
     +                      ABY,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +


             CENSQY (1)=3
             CENSQY (2)=1


             DO N = 1,NDERY

                NINT = MGIJCEN * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DY (M)
                END DO

                CENTER = CENSQY (N)
c                       IF(M.eq.1) CENTER=3
c                       IF(M.eq.2) CENTER=1
         

                IF (CENTER.EQ.1) THEN

                    DEROPA = 1
                    DEROPB = IXAEQB
                    DEROPC = 0
                    IF (CEQA) THEN
                        DEROPC = 1
                        IXC = IXCEQA
                        CY = YN (NUCCEN (IXCEQA))
                    END IF

                ELSE IF (CENTER.EQ.2) THEN

                    DEROPA = IXAEQB
                    DEROPB = 1
                    DEROPC = 0
                    IF (CEQB) THEN
                        DEROPC = 1
                        IXC = IXCEQB
                        CY = YN (NUCCEN (IXCEQB))
                    END IF

                ELSE

                    DEROPA = 0
                    DEROPB = 0
                    DEROPC = 1
                    IXC = 1
                    CY = YN (NUCCEN (1))

                END IF

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB


                IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                    WRITE (*,*) ' Problems forming y-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__nai_derv_pcgto_block '
                    STOP
                END IF


                CALL    OED__NAI_1D_DERV_INTEGRALS
     +
     +                       ( MIJ,MGQPIJ,MGIJCEN,
     +                         NGQP,NGQPCEN,
     +                         SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,DEROPC,
     +                         EXP2A,EXP2B,EXP2AB,
     +                         IXC,PY,CY,
     +                         RTS,
     +                         R1Y,
     +                         BATCH,
     +
     +                                  INT1DY )
     +
     +
             END DO

         ELSE
            
             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Y,R2,
     +                      ABY,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +
         END IF
C---------------------------------------------------------------
C                                                              |
C             ...assemble the 1DZ AB integrals and perform     |
C                differentiation sequence, if necessary.       |
C                                                              |
C---------------------------------------------------------------
         IF (DIFFZ) THEN

             DERA = DERAZ
             DERB = DERBZ

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Z,R2,
     +                      ABZ,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DZ )
     +
     +

             CENSQZ (1)=3
             CENSQZ (2)=1


             DO N = 1,NDERZ

                NINT = MGIJCEN * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DZ (M)
                END DO

                CENTER = CENSQZ (N)
c                  IF(M.eq.1) CENTER=3
c                  IF(M.eq.2) CENTER=1



                IF (CENTER.EQ.1) THEN

                    DEROPA = 1
                    DEROPB = IXAEQB
                    DEROPC = 0
                    IF (CEQA) THEN
                        DEROPC = 1
                        IXC = IXCEQA
                        CZ = ZN (NUCCEN (IXCEQA))
                    END IF

                ELSE IF (CENTER.EQ.2) THEN

                    DEROPA = IXAEQB
                    DEROPB = 1
                    DEROPC = 0
                    IF (CEQB) THEN
                        DEROPC = 1
                        IXC = IXCEQB
                        CZ = ZN (NUCCEN (IXCEQB))
                    END IF

                ELSE

                    DEROPA = 0
                    DEROPB = 0
                    DEROPC = 1
                    IXC = 1
                    CZ = ZN (NUCCEN (1))

                END IF

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB


c--------------------------------
c--------------------------------


                IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                    WRITE (*,*) ' Problems forming z-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__nai_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__NAI_1D_DERV_INTEGRALS
     +
     +                       ( MIJ,MGQPIJ,MGIJCEN,
     +                         NGQP,NGQPCEN,
     +                         SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,DEROPC,
     +                         EXP2A,EXP2B,EXP2AB,
     +                         IXC,PZ,CZ,
     +                         RTS,
     +                         R1Z,
     +                         BATCH,
     +
     +                                  INT1DZ )
     +
     +
             END DO

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Z,R2,
     +                      ABZ,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DZ )
     +
     +
         END IF
C---------------------------------------------------------------------
C                                                                    |
C             ...assemble the 1D AB integrals to the [A|B] batch.    |
C                                                                    |
C---------------------------------------------------------------------
         IF (SHELLP.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_00_F12
     +
     +                    ( MIJ,NGQPCEN,MGIJCEN,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLB.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_A0_F12
     +
     +                    ( SHELLA,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZA,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,R1Y,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLA.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_A0_F12
     +
     +                    ( SHELLB,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,R1Y,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
     +
         ELSE

             CALL    OED__NAI_DERV_INT1D_TO_AB_F12
     +
     +                    ( SHELLA,SHELLB,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZA,NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,R1Y,
     +                      SCALE,RTS,UVAL,STGEXP,
     +
     +                                BATCH )
     +
     +
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
