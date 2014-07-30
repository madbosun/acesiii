      SUBROUTINE A2EVAL_INTS(NTOTPRIM, NTOTCRF, IANGTYPE, 
     &                       XPOINT, YPOINT, ZPOINT, EXPS, PCOEF, 
     &                       PRDTINT, CNT_COORD, TMP1, TMP2, TMP3)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION CNTMU(3), CNTNU(3), EXPS(NTOTPRIM), 
     &          PRDTINT(NTOTCRF, NTOTCRF), PCOEF(NTOTPRIM, NTOTCRF),
     &          CENTER(3),TMP1(NTOTPRIM*NTOTPRIM),
     &          TMP2(NTOTPRIM*NTOTPRIM),TMP3(NTOTPRIM*NTOTPRIM),
     &          IANGTYPE(NTOTPRIM), CNT_COORD(3, NTOTPRIM)
C
      COMMON /HIGHL/ LMNVAL(3, 84), ANORM(84)
      DATA   DZERO /0.0D+00/, ONE /1.00D0/
C
C Loop over number of primitive functions
C
      INDX = 0
C
      CENTER(1) = XPOINT
      CENTER(2) = YPOINT
      CENTER(3) = ZPOINT
C
      DO LFTPRIM = 1, NTOTPRIM
C
         INDX  = INDX + 1 
         ITYPE = IANGTYPE(INDX)
         DO IXYZ = 1, 3
            CNTMU(IXYZ) = CNT_COORD(IXYZ, INDX)
         ENDDO 
C                     
         JNDX = 0
         DO RGTPRIM = 1, NTOTPRIM
C                  
            JNDX  = JNDX + 1 
            JTYPE = IANGTYPE(JNDX)
            DO IXYZ = 1, 3
               CNTNU(IXYZ) = CNT_COORD(IXYZ, JNDX)
            ENDDO 
C
C We can built the product here!
C
                  CALL A2BULT_PRDUCT(INDX, JNDX, ITYPE, JTYPE, 
     &                               NTOTPRIM, EXPS(INDX), 
     &                               EXPS(JNDX), CENTER,
     &                               CNTMU, CNTNU, TMP2)
C
C Loop over primitives end here!
C
         ENDDO
      ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) "The product integral in primitive basis"
      CALL OUTPUT(TMP2, 1, NTOTPRIM, 1, NTOTPRIM, NTOTPRIM, 
     &            NTOTPRIM, 1)
#endif
C
#ifdef _DEBUG_LVLM2
       Write(6,*)
       Write(6,*) "@-a2eval_ints CONTRACTION COEFICIENTS"
       CALL OUTPUT(PCOEF, 1, NTOTPRIM, 1, NTOTCRF, NTOTPRIM,
     &                  NTOTCRF, 1)
#endif
C 
C Built the contracted product functions for this shell.
C
      CALL XGEMM('N', 'N', NTOTPRIM, NTOTCRF, NTOTPRIM, ONE, TMP2,
     &            NTOTPRIM, PCOEF, NTOTPRIM, DZERO, TMP3, NTOTPRIM)
C
      CALL XGEMM('T', 'N', NTOTCRF, NTOTCRF, NTOTPRIM, ONE, PCOEF, 
     &            NTOTPRIM, TMP3, NTOTPRIM, DZERO, PRDTINT, NTOTCRF)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The product integral in contracted basis"
      CALL OUTPUT(PRDTINT, 1, NTOTCRF, 1, NTOTCRF, NTOTCRF, NTOTCRF, 1)
#endif
C
      RETURN 
      END
