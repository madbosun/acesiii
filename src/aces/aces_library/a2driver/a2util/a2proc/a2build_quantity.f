      SUBROUTINE A2BUILD_QUANTITY(WORK, MAXCOR, PRDUTINT, QUANTITY, 
     &                            SPIN_D, DENSITY_TYPE, 
     &                            NMBR_OF_PERTS, IPICK_PERT, NBFNS,
     &                            NAOBFNS, ISCF_TDEN, ICOR_TDEN,  
     &                            ISCF_DDEN, ICOR_DDEN, 
     &                            IBEGIN_P_DENS, MAX_GRID_POINTS, 
     &                            POST_SCF, GRID_TYPE, KCUBE, IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPIN_D, POST_SCF
C
      CHARACTER*7 DENSITY_TYPE
      CHARACTER*2 GRID_TYPE
C
      DIMENSION WORK(MAXCOR), PRDUTINT(NAOBFNS, NAOBFNS),
     &          QUANTITY(0:4, MAX_GRID_POINTS)
      DATA TWO /2.0D0/
C
      IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
         SCF_TOT_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                          WORK(ISCF_TDEN), 1)
         COR_TOT_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                         WORK(ICOR_TDEN), 1)
C
         IF (.NOT. SPIN_D) THEN
C
CSSS            DO IND_CUBE = 1, KCUBE
               QUANTITY(0, KCUBE) = SCF_TOT_DENSITY
               IF (POST_SCF)  QUANTITY(1, KCUBE) = 
     &                                COR_TOT_DENSITY
CSSS            END DO
         ENDIF
C
         IF (IUHF .NE. 0 .AND. SPIN_D) THEN
            SCF_DIF_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT,
     &                             1, WORK(ISCF_DDEN), 1)
            COR_DIF_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT,
     &                             1, WORK(ICOR_DDEN), 1)
C
CSSS            DO IND_CUBE = 1, KCUBE
               QUANTITY(0, KCUBE) = SCF_DIF_DENSITY
               IF (POST_SCF)  QUANTITY(1, KCUBE) = 
     &                                 COR_DIFF_DENSITY
CSSS            END DO
         ENDIF
C
      ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN

          DO IPERT = 1, NMBR_OF_PERTS
             IF (IPERT .EQ. IPICK_PERT) THEN
            
                I_TOT_PDEN_N = IBEGIN_P_DENS + 2*(IPERT-1)*
     &                         NAOBFNS*NAOBFNS
                I_DIF_PDEN_N =  I_TOT_PDEN_N + NAOBFNS*NAOBFNS

                TOT_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                            WORK(I_TOT_PDEN_N), 1)
                DIF_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                            WORK(I_DIF_PDEN_N), 1)
C
CSSS                DO IND_CUBE = 1, KCUBE
                   QUANTITY(0, KCUBE) = TOT_PDEN_N*TWO
                   QUANTITY(1, KCUBE) = DIF_PDEN_N*TWO
CSSS                END DO
C
             END IF
          END DO
c
      ELSE IF (DENSITY_TYPE .EQ. "DEFINED") THEN

         OPEN (UNIT=30, FILE='VPOUT', FORM='UNFORMATTED', 
     &         STATUS='OLD')
         IRWND = 0
         IERR  = 0
         IOFF  = 1
C
         DO IPERT = 1, NMBR_OF_PERTS

            CALL SEEKLB ('   DEN  ', IERR, IRWND, 30)
            IF (IERR .NE. 0) CALL ERREX
C
            CALL LOAD_DEFINED(WORK(IOFF), NAOBFNS*NAOBFNS, NAOBFNS)
C
C#ifdef _DEBUG_LVLM2
            CALL HEADER ('Matrix elements of the defined operator', 1
     &                   , 6)
            CALL OUTPUT(WORK(IOFF), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &                  NAOBFNS, 1)
C#endif 
            IRWND = 1
            IOFF  = IOFF + NAOBFNS*NAOBFNS
            
         END DO
C
         DO IPERT = 1, NMBR_OF_PERTS
             IF (IPERT .EQ. IPICK_PERT) THEN

                I_TOT_PDEN_N = IBEGIN_P_DENS + (IPERT-1)*
     &                         NAOBFNS*NAOBFNS

                TOT_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                            WORK(I_TOT_PDEN_N), 1)
C
CSSS                DO IND_CUBE = 1, KCUBE
                   QUANTITY(0, KCUBE) = TOT_PDEN_N*TWO
CSSS                END DO
C
             END IF
          END DO
c
      Print*, "Quantity of interest (perturbed total and/or spin",
     &        " density)"
      Write(6,*), KCUBE
CSS      Do IND_CUBE =1, KCUBE
      Write(*,'((1x,2F10.6))') QUANTITY(0,KCUBE),
     &                         QUANTITY(1,KCUBE)
CSSS      Enddo
      ENDIF
C
      RETURN
      END 
C  
