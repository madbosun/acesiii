       CHARACTER*(*) FUNCTION FILTSP (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Filter spaces - reduce each run of whitespace to a 
C               single space character.
C
C Arguments:    STRING      character string (input only)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C $Log: filtsp.f,v $
C Revision 1.1.1.1  2003/04/02 19:21:35  aces
C INITIAL 2.4 IMPORT
C
c Revision 4.1  89/08/19  14:06:48  bernhold
c From the strings library, but this one works right!
c 
C    Revision 1.1  88/01/14  13:14:25  bernhold
C    Initial revision
C    
C
C System:       Standard FORTRAN 77
C 
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       CHARACTER*(*) STRING
       INTEGER I, J
       CHARACTER*1 C
       LOGICAL INSPAC
C
C      I is position in input string, J is position in output string
C
       J = 1
       INSPAC = .FALSE.
       DO 100 I = 1, LEN(STRING)       
          C = STRING(I:I)
C
C         Whitespace is defined as in the routine ISSPAC, to be:
C            HT, LF, VT, FF, CR, or space
          IF ((LGE(C,CHAR(9)) .AND. LLE(C,CHAR(13))) .OR. C.EQ.' ')
     1    THEN
C
C            Copy only the first space in a run of them
             IF (INSPAC) THEN
                FILTSP(J:J) = ' '
                INSPAC = .TRUE.
                J = J + 1
             ENDIF
C
C         This isn't a space, so just copy it.
          ELSE
             INSPAC = .FALSE.
             FILTSP(J:J) = C
             J = J + 1
          ENDIF
  100  CONTINUE
C
       RETURN
       END
