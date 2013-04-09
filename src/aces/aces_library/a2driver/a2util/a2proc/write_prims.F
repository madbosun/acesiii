       Subroutine Write_Prims(NATOMS, NTOT_PRIMS, NBFNS, NAOBFNS,
     &                        COORD, ALPHA, PCOEF, NUNQSHL, NSHL,
     &                        NANGMOMSHL, NCONFUNSHL, NPRIMFUNSHL,
     &                        CNT_COORD, ANGTYPE, IATMCHRG)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "flags.h"
#include "machsp.com"
#include "jodaflags.com"                   
#include "mxatms.par"
C#include "istart.com"                      
#include "baslims.par"
C    
      DIMENSION ALPHA(NTOT_PRIMS), IPOPF(NATOMS), COORD(NATOMS*3), 
     &          NSHL(NATOMS), NANGMOMSHL(NATOMS, NUNQSHL),
     &          NCONFUNSHL(NATOMS, NUNQSHL), 
     &          NPRIMFUNSHL(NATOMS,NUNQSHL), PCOEF(NTOT_PRIMS, NBFNS), 
     &          IANGTYPE(NTOT_PRIMS), ATOM_XYZ(4, MXATMS),
     &          CNT_COORD(3,NTOT_PRIMS),CNTMU(3)
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
C
      DATA XTANG /0.529177249D0/ 
C
#ifdef _DEBUG_LVL0
       Write(6,*) 
       Write(6,"(a,5I4)"), "Input varaibles:",  Natoms, 
     &         nunqshl, nbfns, NTOT_PRIMS 
       Write(6,*)
       write(*, '(a)'), "The exponents"
       Write(*, '(4(1x,F10.5))') (alpha(i), i=1, NTOT_PRIMS)
       Write(6,*)
       Write(*, '(a)'), "Cartesian Coords."
       Write(*, '(3F10.5)') (Coord(i), i=1, 3*NATOMS)
       Write(6,*)
       Write(*, '(a)'), "The NSHL array"
       Write(*, '(4I5)') (Nshl(i), i=1, NATOMS)
       Write(6,*)
       Write(6,*) "The number of primitive in a shell"
       Write(6, '(4I5)') ((NPRIMFUNSHL(i,j), J=1, Nshl(i)),  
     &                    i=1,NATOMS)
       Write(6,*)
       Write(6,*) "The number of contracted functions in a shell"
       Write(6, '(4I5)') ((NCONFUNSHL(i,j), J=1, Nshl(i)),  
     &                    i=1,NATOMS)
       Write(6,*)
       Write(6,*) "The angular momentum of shells"
       Write(6, '(4I5)') ((NANGMOMSHL(i,j), J=1, Nshl(i)), 
     &                    I=1,NATOMS)
       write(6,*)
       Write(*, '(a)') "The Contractions Coef."
       write(*, '(6F10.5)') ((pcoef(i, j), i=1, NTOT_PRIMS), 
     &                        j=1,nbfns)
#endif              
      OPEN(UNIT=IUNIT, FILE="BASIS_LOG", FORM="FORMATTED", STATUS="NEW")
C       
      CALL A2GET_ANG_TYPE(NATOMs, NTOT_PRIMS, NUNQSHL, NSHL, NANGMOMSHL, 
     &                    NPRIMFUNSHL, IANGTYPE, COORD, CNT_COORD)
C
      INDEX = 0
      DO LFTPRIM = 1, NTOTPRIM
C
         INDX  = INDX + 1
         ITYPE = IANGTYPE(INDX)
C
         DO IXYZ = 1, 3
            CNTMU(IXYZ) = CNT_COORD(IXYZ, INDX)
         ENDDO
C
         LI = LMNVAL(1, ITYPE)
         MI = LMNVAL(2, ITYPE)
         NI = LMNVAL(3, ITYPE)

      ENDDO 
C

      Return
      END
		
		
	
