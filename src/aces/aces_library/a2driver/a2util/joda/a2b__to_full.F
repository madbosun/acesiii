         SUBROUTINE  A2B__TO_FULL
     +
     +                    ( NATOMS,
     +                      IMAX,  ZMAX,
     +                      ICORE, ZCORE )
     +
     +
C-----------------------------------------------------------------------
C  OPERATION   : A2B__TO_FULL
C  MODULE      : ACESII BASIS CONVERTER
C  MODULE-ID   : A2B
C  SUBROUTINES : A2B__FNC_SIZE
C                A2B__READ_REST
C                A2B__WRITE_NEW
C  DESCRIPTION : This routine will read the ZMAT.BAS file made by JODA
C                and generate a new formatted basis set file where each
C                contracted function is written individually in all its
C                glory!
C
C                Input:
C
C                     NATOMS = The number of atoms in the ZMAT.BAS
C                     IMAX   = The amount of integer memory
C                     ZMAX   = The amount of floating point memory
C                     ICORE  = The array of memory for integers
C                     ZCORE  = The array of memory for floating point
C                              numbers
C
C                Output:
C
C                     There is NO output for this routine in terms of
C                     variables!  A new file is generated though, called
C
C                                      ZMAT.NEW
C
C
C  AUTHOR      : Thomas Watson Jr.
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    STRLN,JUNK
         INTEGER    I,J,IATOM
         INTEGER    UACES,UNEW
         INTEGER    IMAX,ZMAX
         INTEGER    NATOMS,NSHELL
         INTEGER    MAXCNTR,MAXPRIM

         LOGICAL    LOPEN, LEXIST, CLSE

         CHARACTER*80  TITLE1,TITLE2

         INTEGER           ICORE (1:IMAX)
         DOUBLE PRECISION  ZCORE (1:ZMAX)

         INTEGER    ILANG,ICNTR,IPRIM,ISCR,IEND
         INTEGER    ZSCR1,ZSCR2,ZEXPN,ZCOEF,ZEND

         PARAMETER  ( UACES = 12 )
         PARAMETER  ( UNEW  = 46 )
C
C
C------------------------------------------------------------------------
C
C
C             ...Carefully check and for existence and open ZMAT.BAS
C
C
         INQUIRE (FILE = 'GENB.TMP', OPENED = LOPEN, EXIST = LEXIST )

         IF (.NOT. LEXIST) THEN
            WRITE (*,*) ' @A2B__TO_FULL: GENB.TMP does not exist!'
            WRITE (*,*) '                Cannot format a file that ',
     +                                   'does not exist!'
         ENDIF

         IF (.NOT. LOPEN) THEN
            OPEN (UACES, FILE = 'GENB.TMP', STATUS = 'OLD')
            CLSE = .TRUE.
            REWIND (UACES)
         ELSE
            REWIND (UACES)
            CLSE = .FALSE.
         ENDIF
C
C
C             ...Carefully check and for existence and open the new
C                formatted ZMAT.NEW
C
C
         INQUIRE (FILE = 'ZMAT.BAS', EXIST = LEXIST )

         IF (LEXIST) THEN
            WRITE (*,*) ' @A2B__TO_FULL: The newly formatted ZMAT.BAS'
            WRITE (*,*) '                already exists!'
            RETURN
         ELSE
            OPEN   (UNEW, FILE = 'ZMAT.BAS', STATUS = 'NEW')
            REWIND (UNEW)
         ENDIF

C
C
C             ...Loop over the number of atoms in the ZMAT.BAS
C
C
         DO IATOM = 1, NATOMS
C
C
C             ...Read in the ZMAT.BAS information.  We need it all
C                to properly format the file!
C
C
            READ (UACES,'(A)')  TITLE1     ! Element symbol and basis set
            READ (UACES,'(A)')  TITLE2     ! Where basis set is found
            READ (UACES,'(A)')  JUNK       ! Blank line
            READ (UACES,'(I3)') NSHELL     ! Number of shells
            STRLN = LEN_TRIM (TITLE1)      ! Trim the basis set name only
C
C
C             ...Equipped with NSHELL, we need to call a subroutine
C                to continue reading the next three lines!
C
C
            ILANG = 1
            ICNTR = ILANG + NSHELL
            IPRIM = ICNTR + NSHELL
            IEND  = IPRIM + NSHELL

            IF ((IEND - ILANG) .GT. IMAX) THEN
               WRITE (*,*) ' @A2B__TO_FULL: Insufficient integer ',
     +                                     'memory.  No output will ',
     +                                     'be produced.'
               RETURN
            END IF
 
            CALL  A2B__FNC_SIZE
     +
     +                 ( UACES,NSHELL,
     +                   MAXCNTR,MAXPRIM,
     +
     +                           ICORE (ILANG),
     +                           ICORE (ICNTR),ICORE (IPRIM) )
     +
     +
C
C
C             ...Now that we have the angular momentum, the number
C                of contracted functions, and the number of primitive
C                functions, we can read the rest of the ZMAT.BAS
C                for this atom.
C
C
            ZEXPN = 1
            ZCOEF = ZEXPN + NSHELL * MAXPRIM
            ZSCR1 = ZCOEF + NSHELL * MAXCNTR * MAXPRIM
            ZSCR2 = ZSCR1 + NSHELL * MAXCNTR * MAXPRIM
            ZEND  = ZSCR2 + NSHELL * MAXCNTR * MAXPRIM

            IF ((ZEND - ZEXPN) .GT. ZMAX) THEN
               WRITE (*,*) ' @A2B__TO_FULL: Insufficient fl. pt. ',
     +                                     'memory.  No output will ',
     +                                     'be produced.'
               RETURN
            END IF

            CALL  A2B__READ_REST
     +
     +                 ( UACES,
     +                   NSHELL, MAXCNTR, MAXPRIM,
     +                   ICORE (ILANG),ICORE (ICNTR),
     +                   ICORE (IPRIM),ZCORE (ZSCR1),
     +
     +                                       ZCORE (ZEXPN),
     +                                       ZCORE (ZCOEF) )
     +
     +
C
C
C             ...Now we have the basis set information for a
C                given atom.  We can proceed to write out the
C                newly formated basis set in ZMAT.NEW!
C
C
            ISCR = IPRIM + NSHELL
            IEND = ISCR  + NSHELL * MAXCNTR * MAXPRIM

            IF ((IEND - ILANG) .GT. IMAX) THEN
               WRITE (*,*) ' @A2B__TO_FULL: Insufficient integer ',
     +                                     'memory.  No output will ',
     +                                     'be produced.'
               RETURN
            END IF


            CALL  A2B__WRITE_NEW
     +
     +                 ( UNEW,TITLE1(1:STRLN),TITLE2,
     +                   NSHELL,MAXCNTR,MAXPRIM,
     +                   ICORE (ILANG),
     +                   ICORE (ICNTR),ICORE (IPRIM),
     +                   ICORE (ISCR) ,ZCORE (ZSCR1),
     +                   ZCORE (ZSCR2),
     +                   ZCORE (ZEXPN),ZCORE (ZCOEF) )
     +
     +
C
C
C             ...End the loop going over the number of atoms!
C
C
         END DO
C
C
C             ...If the ZMAT.BAS file was closed upon entering,
C                close it now!  Also close the ZMAT.NEW file!
C
C
         IF (CLSE) CLOSE (UACES, STATUS = 'DELETE')
         CLOSE (UNEW, STATUS = 'KEEP')
C
C
C             ...ready!
C
C
         RETURN
         END

