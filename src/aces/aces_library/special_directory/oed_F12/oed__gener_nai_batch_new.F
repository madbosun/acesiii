         SUBROUTINE  OED__GENER_NAI_BATCH_NEW
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      NUCLEI,
     +                      XN,YN,ZN,NCHARGE,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      SPHERIC,
     +                      SCREEN,
     +                      ICORE,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__GENER_NAI_BATCH
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__NAI_CSGTO
C  DESCRIPTION : Main operation that drives the calculation of a batch
C                of contracted electron nuclear attraction integrals.
C
C
C                  Input (x = 1 and 2):
C
C                    IMAX,ZMAX    =  maximum integer + flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                    NUCLEI       =  # of nuclear attraction centers
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    NCHARGE      =  the nuclear charges for all
C                                    nuclear attraction centers
C                    ALPHA        =  primitive exponents for csh
C                                    1 and 2 in that order
C                    CC           =  contraction coefficient for csh
C                                    1 and 2 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1 and 2 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C                    SCREEN       =  is true, if screening will be
C                                    done at primitive integral level
C                    ICORE        =  integer scratch space
C                    ZCORE (part) =  flp scratch space
C
C
C                  Output:
C
C                    NBATCH       =  # of integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    integral
C                    ZCORE        =  full batch of contracted (1|2)
C                                    nuclear attraction integrals over
C                                    cartesian or spherical gaussians
C                                    starting at ZCORE (NFIRST)
C
C
C
C                            !!! IMPORTANT !!!
C
C                For performance tuning, please see the include file
C                'oed__tuning.inc'.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C                      >>>  COPYRIGHT NOTICE  <<<
C                             Program OED
C
C    The program package OED is copyrighted 2004 by N. Flocke.
C    All rights are reserved, except as noted herein. OED was produced
C    in part under the auspices of the U. S. Department of Defense.
C    The U. S. Government retains a nonexclusive, royalty-free license
C    to publish or reproduce any published form of OED, or to allow
C    others to do so, for U. S. Government purposes. This copy of OED
C    may not be redistributed nor reproduced, in whole, in part, or as
C    part of another code package without the prior, explicit, written
C    permission from the author. This copy of OED may be used for
C    scientific research provided that any and all published papers
C    and/or presentations of results generated using OED shall
C    acknowledge that fact and the authorship of the code specifically
C    and accurately.
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INCLUDE     'oed__ftable.inc'
         INCLUDE     'oed__tuning.inc'

         LOGICAL     SCREEN
         LOGICAL     SPHERIC

         INTEGER     IMAX,ZMAX
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NUCLEI
         INTEGER     SHELL1,SHELL2

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2

         DOUBLE PRECISION  XN      (1:NUCLEI)
         DOUBLE PRECISION  YN      (1:NUCLEI)
         DOUBLE PRECISION  ZN      (1:NUCLEI)
         DOUBLE PRECISION  NCHARGE (1:NUCLEI)

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)
C
C
C------------------------------------------------------------------------
C
C
C             ...call csgto routine.
C
C
         CALL  OED__NAI_CSGTO
     +
     +              ( IMAX,ZMAX,
     +                NALPHA,NCOEFF,NCSUM,
     +                NCGTO1,NCGTO2,
     +                NPGTO1,NPGTO2,
     +                SHELL1,SHELL2,
     +                X1,Y1,Z1,X2,Y2,Z2,
     +                NUCLEI,
     +                XN,YN,ZN,NCHARGE,
     +                ALPHA,CC,CCBEG,CCEND,
     +                FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                L1CACHE,TILE,NCTROW,
     +                SPHERIC,SCREEN,
     +                ICORE,
     +
     +                          NBATCH,
     +                          NFIRST,
     +                          ZCORE )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
