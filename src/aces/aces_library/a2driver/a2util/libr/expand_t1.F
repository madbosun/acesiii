
      SUBROUTINE EXPAND_T1(T1_EXPND, T1_DROP, NIRREP, NOCC_EXPND,
     &                     NVRT_EXPND, NOCC, NVRT, NOCCO_EXPND,
     &                     NVRTO_EXPND, SYMP_OVFUL)
C
C This will expand the symmetry-packed T1 space from the drop-mo to full-mo by
C adding zeros to the excitations that involve dropped mos: t(i,a)=0
C for all inactive i and a. Ajith Perera, 07/2005.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER SYMP_OVFUL
      DIMENSION T1_EXPND(NVRTO_EXPND*NOCCO_EXPND), T1_DROP(NVRT*NOCC),
     &          NOCC_EXPND(8), NVRT_EXPND(8)
      COMMON /SYMDROP/ NDRPOP(8),NDRVRT(8)

      IOFF_FULL = 0
      IOFF_ACT  = 0

      CALL ZERO(T1_EXPND, NVRTO_EXPND*NOCCO_EXPND)
      CALL GETREC(20,'JOBARC','NDROPPOP',NIRREP,NDRPOP)
      CALL GETREC(20,'JOBARC','NDROPVRT',NIRREP,NDRVRT)

      DO IRREP = 1, NIRREP
             NVRT = NVRT_EXPND(IRREP)
             NOCC = NOCC_EXPND(IRREP)
         NVRT_DRP = NDRVRT(IRREP)
         NOCC_DRP = NDRPOP(IRREP)
         NVRT_ACT = NVRT - NVRT_DRP
         NOCC_ACT = NOCC - NOCC_DRP

         IOFF_FULL = IOFF_FULL + NVRT*NOCC
         IOFF_ACT  = IOFF_ACT  + NVRT_ACT*NOCC_ACT
      END DO
      SYMP_OVFUL = IOFF_FULL
      IOFF_FULL  = IOFF_FULL + 1
      IOFF_ACT   = IOFF_ACT  + 1

      DO IRREP = NIRREP, 1, -1
             NVRT = NVRT_EXPND(IRREP)
             NOCC = NOCC_EXPND(IRREP)
         NVRT_DRP = NDRVRT(IRREP)
         NOCC_DRP = NDRPOP(IRREP)
         NVRT_ACT = NVRT - NVRT_DRP
         NOCC_ACT = NOCC - NOCC_DRP

         IOFF_FULL = IOFF_FULL - NVRT*NOCC
         IOFF_ACT  = IOFF_ACT  - NVRT_ACT*NOCC_ACT

         CALL FULL_T1(T1_EXPND(IOFF_FULL), T1_DROP(IOFF_ACT),
     &                NVRT_ACT, NOCC_ACT, NVRT, NOCC,
     &                NVRT_DRP, NOCC_DRP)

      END DO

      RETURN
      END

