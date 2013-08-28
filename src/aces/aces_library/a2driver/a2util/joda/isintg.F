       LOGICAL FUNCTION ISINTG (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      True if STRING is a valid integer number
C
C Arguments:    STRING   character string (input only)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C $Log: isintg.f,v $
C Revision 1.1.1.1  2003/04/02 19:21:35  aces
C INITIAL 2.4 IMPORT
C
c Revision 4.1  89/08/19  14:07:30  bernhold
c From the strings library, but this one works right!
c 
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       CHARACTER*(*) STRING
       CHARACTER*1 C
       INTEGER I, J, linblnk
C
C      Cases: No string => not integer (or anything else)
C             1 char.   => must be a digit to be an integer
C             >1 char.  => first char can be +/- or digit
C
       IF (LEN(STRING) .EQ. 0) THEN
          ISINTG = .FALSE.
       ELSEIF (LEN(STRING) .EQ. 1) THEN
          ISINTG= .TRUE.
          IF (LLT(STRING,'0') .OR. LGT(STRING,'9')) ISINTG = .FALSE.
       ELSE
C            First, find the last non-blank character - an linblnk without
C            actually calling the routine (keep these basic routines indep.)
          DO 10 I = LEN(STRING), 1, -1
             IF (STRING(I:I) .NE. ' ') THEN
                linblnk = I
                GOTO 20
             ENDIF
 10       CONTINUE
c     
c  string is blank -- not a number
          isintg = .false.
          return
c
C            I is used to start the loop for digits at the proper pos'n
 20       I = 1
 100      IF (STRING(I:I) .EQ. ' ') THEN
             I = I + 1
             GOTO 100
          ENDIF
C
          IF (STRING(I:I) .EQ. '-' .OR. STRING(I:I) .EQ. '+') THEN
             I = I + 1
          ENDIF
C            All remaining places must be digits
          ISINTG = .TRUE.
          DO 200 J = I, linblnk
             C = STRING(J:J)
             IF (LLT(C,'0') .OR. LGT(C,'9')) ISINTG = .FALSE.
 200      CONTINUE
       ENDIF
       RETURN
       END
