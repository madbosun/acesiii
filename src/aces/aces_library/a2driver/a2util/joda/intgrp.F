      INTEGER FUNCTION INTGRP(PGRP)
C
C   GET THE INTERNAL NUMBER OF THE IRREP
C  
      CHARACTER*3 IGROUP,PGRP
      DIMENSION IGROUP(8)
      DATA IGROUP /'C1 ','C s','C2 ','C i','D2 ','C2v','C2h','D2h'/
C
      DO 1 I=1,8
         IF(IGROUP(I).EQ.PGRP) GOTO 2
 1    CONTINUE
      WRITE(6,3) 
 3    FORMAT(T3,'@INTGRP-F Symbol of point group do not match')
      CALL ERREX
 2    INTGRP=I
      RETURN
      END
