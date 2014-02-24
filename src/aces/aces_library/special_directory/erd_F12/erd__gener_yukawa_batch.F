         SUBROUTINE  ERD__GENER_YUKAWA_BATCH
     +
     +                    ( IMAX,ZMAX,STGEXP,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
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
C  OPERATION   : ERD__GENER_ERI_BATCH
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__1111_CSGTO
C                ERD__CSGTO
C  DESCRIPTION : Main operation that drives the calculation of a batch
C                of contracted electron repulsion integrals.
C
C
C                  Input (x = 1,2,3 and 4):
C
C                    IMAX,ZMAX    =  maximum integer,flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2,3 and 4
C                    ALPHA        =  primitive exponents for csh
C                                    1,2,3,4 in that order
C                    CC           =  contraction coefficient for csh
C                                    1,2,3,4 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1,2,3,4 in that order. They are
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
C                    ZCORE        =  full batch of contracted (12|34)
C                                    integrals over cartesian or
C                                    spherical gaussians starting at
C                                    ZCORE (NFIRST)
C
C
C
C                            !!! IMPORTANT !!!
C
C                For performance tuning, please see the include file
C                'erd__tuning.inc'.
C                
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INCLUDE     'erd__ftable.inc'
         INCLUDE     'erd__tuning.inc'

         LOGICAL     SCREEN
         LOGICAL     SPHERIC

         INTEGER     IMAX,ZMAX,ii
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         double precision STGEXP

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)

c           do ii=1,NCOEFF
c           write(6,*) CC(ii)
c           enddo
C
C------------------------------------------------------------------------
C
C
C
C
C
C---------------Generation of the batch of Yukawa integrals--------------

             CALL  ERD__CSGTO_YUKAWA
     +
     +                  ( IMAX,ZMAX,STGEXP,
     +                    NALPHA,NCOEFF,NCSUM,
     +                    NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                    NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                    SHELL1,SHELL2,SHELL3,SHELL4,
     +                    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                    ALPHA,CC,CCBEG,CCEND,
     +                    FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                    L1CACHE,TILE,NCTROW,
     +                    SPHERIC,SCREEN,
     +                    ICORE,
     +
     +                              NBATCH,
     +                              NFIRST,
     +                              ZCORE )
     +
     +






C             ...ready!
C
C
         RETURN
         END
