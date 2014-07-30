      SUBROUTINE A2GET_ANG_TYPE(MAXATMS, NTOTPRIM, MAXSHELL, NSHL, 
     &                          NANGMOMSHL, NPRIMFUNSHL, IANGTYPE,
     &                          COORD, CNT_COORD)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION NANGMOMSHL(MAXATMS, MAXSHELL), IP(20), FAC(9,9), 
     &          IOFFST(15), IANGTYPE(NTOTPRIM), NSHL(MAXATMS),
     &          NPRIMFUNSHL(MAXATMS,MAXSHELL), 
     &          CNT_COORD(3, NTOTPRIM), COORD(3*MAXATMS)
C      
      CALL SETRHF(FAC, IOFFST, IP)
C
      NFUNC_COUNT = 0
      DO IATMS = 1, MAXATMS
C                  
            DO ISHL = 1, NSHL(IATMS)
C
               IANGMOM = NANGMOMSHL(IATMS, ISHL)
               IPRIM   = NPRIMFUNSHL(IATMS, ISHL)
               IINTTYP = IOFFST(IANGMOM)
               ITYPE   = IINTTYP
C
               DO IMOMN = 1, IANGMOM
                  ITYPE = ITYPE + 1
                  DO NPRIM = 1, IPRIM
                     NFUNC_COUNT = NFUNC_COUNT + 1
                     IANGTYPE(NFUNC_COUNT) = ITYPE
C
                     ICOORD = (IATMS - 1)*3
                     DO IXYZ = 1, 3
                        CNT_COORD(IXYZ, NFUNC_COUNT) = 
     &                                  COORD(ICOORD + IXYZ)
                     ENDDO
   
                  ENDDO
               ENDDO
C
            ENDDO
      ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,I4)") "The number of primitive functions", 
     &                   NFUNC_COUNT 
      Write(6,*)
      write(6,"(a)") "The angular momentum type of each function"
      Write(6,"(6i4)") (IANGTYPE(i), i=1, NFUNC_COUNT)
#endif
      RETURN
      END
 

