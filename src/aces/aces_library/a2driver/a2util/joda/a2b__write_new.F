         SUBROUTINE  A2B__WRITE_NEW
     +
     +                    ( UNEW,TITLE1,TITLE2,
     +                      NSHELL,MAXCNTR,MAXPRIM,
     +                      LANGSH,NCNTSH,NPRMSH,
     +                      ISCR,MEXP,MCEF,EXPNT,COEFS )
     +
     +
C-----------------------------------------------------------------------
C  OPERATION   : A2B__WRITE_NEW
C  MODULE      : ACESII BASIS CONVERTER
C  MODULE-ID   : A2B
C  SUBROUTINES : none
C  DESCRIPTION : This routine will analyze the given exponents and
C                coefficients, organize and map them, and then generate
C                a new ZMAT.BAS (ZMAT.NEW).
C
C                Input:
C
C                     UACES   = The unit number for the ZMAT.BAS
C                     TITLE1  = The first line of the ZMAT.BAS,
C                               i.e. the atom and basis set name
C                     TITLE2  = The second line of the ZMAT.BAS,
C                               i.e. where the basis set was obtained
C                     NSHELL  = The number of shells in the ZMAT.BAS
C                     MAXCNTR = The maximum number of contracted
C                               functions in the basis set
C                     MAXPRIM = The maximum number of primitive
C                               functions in the basis set
C                     LANGSH  = An array of angular momentum numbers
C                     NCNTSH  = An array of the number of contracted
C                               functions for a given angular momenta
C                     NPRMSH  = An array of the number of primitive
C                               functions for a given angular momenta
C                     ISCR    = An integer scratch array
C                     MEXP    = A modified exponent matrix
C                     MCEF    = A modified coefficient matrix that matches
C                               MEXP
C                     EXPNT   = A matrix of the exponents for a given
C                               angular momenta 
C                     COEFS   = A matrix of the coefficients for a 
C                               given angular momenta and contracted
C
C                Output:
C
C                     There is NO output for this routine in terms of
C                     variables!  A new file is generated though.
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

         INTEGER    I,J,K,L,KK,IANG,JUNK
         INTEGER    ISTART,IEND,PART,LOOP
         INTEGER    IPRIM,ICNTR,ICNT,IOFF,ICOUNT
         INTEGER    UNEW,ACESEXP,ACESCNT,ACESSH
         INTEGER    NSHELL,MAXCNTR,MAXPRIM,TOTSH

         CHARACTER*(*)  TITLE1,TITLE2

         DOUBLE PRECISION  DCHECK
         DOUBLE PRECISION  ZERO,CTHRESH

         INTEGER    LANGSH (1:NSHELL)
         INTEGER    NCNTSH (1:NSHELL)
         INTEGER    NPRMSH (1:NSHELL)

         DOUBLE PRECISION  EXPNT (1:NSHELL,1:MAXPRIM)
         DOUBLE PRECISION  COEFS (1:NSHELL,1:MAXCNTR,1:MAXPRIM)

         INTEGER           ISCR  (1:NSHELL*MAXPRIM*MAXCNTR)
         DOUBLE PRECISION  MEXP  (1:NSHELL,1:MAXCNTR,1:MAXPRIM)
         DOUBLE PRECISION  MCEF  (1:NSHELL,1:MAXCNTR,1:MAXPRIM)

         PARAMETER  ( ZERO    = 0.0D0   )
         PARAMETER  ( CTHRESH = 0.1D-06 )
         PARAMETER  ( ACESEXP = 5       )
         PARAMETER  ( ACESCNT = 7       )
         PARAMETER  ( ACESSH  = 14      )
C
C
C------------------------------------------------------------------------
C
C
C             ...Zero out the appropriate arrays!
C
C
         DO I = 1, NSHELL
         DO J = 1, MAXCNTR
         DO K = 1, MAXPRIM

            MEXP (I,J,K) = ZERO
            MCEF (I,J,K) = ZERO

         END DO
         END DO
         END DO
C
C
C             ...Write out the titles!
C
C
         WRITE (UNEW,'(A)')  TITLE1
         WRITE (UNEW,'(A)')  TITLE2
         WRITE (UNEW,  *  )  ''
C
C
C             ...Write out the newly formatted 'number' of shells and
C                angular momentum line!
C
C
         TOTSH = 0
         IOFF  = 0
         DO I = 1, NSHELL
            ICNT  = NCNTSH (I)
            IANG  = LANGSH (I)
            TOTSH = TOTSH + ICNT
            DO J = 1, ICNT
               ISCR (IOFF + J) = IANG
            END DO
            IOFF = IOFF + ICNT
         END DO

         WRITE (UNEW,'(I3)') TOTSH

         IF (TOTSH .GT. ACESSH) THEN

            PART = MOD (TOTSH,ACESSH)
            LOOP = (TOTSH - PART) / ACESSH
            ISTART = 1
            IEND   = ISTART + ACESSH - 1
            DO I = 1, LOOP
               WRITE (UNEW,'(14I5)') (ISCR (J), J = ISTART, IEND)
               ISTART = ISTART + ACESSH
               IF ((TOTSH - IEND) .GT. ACESSH) THEN
                  IEND = IEND + ACESSH
               ELSE
                  IEND = TOTSH
               ENDIF
            END DO

            IF (PART .GT. 0) THEN
               WRITE (UNEW,'(14I5)') (ISCR (J), J = ISTART, IEND)
            END IF

         ELSE
            WRITE (UNEW,'(14I5)') (ISCR (J), J = 1, TOTSH)
         END IF
C
C
C             ...This is the easy part.  Since each contracted function
C                is written out in all its glory, the number of 
C                contracted functions is 1!
C
C
         IOFF  = 0
         DO I = 1, NSHELL
            ICNT  = NCNTSH (I)
            DO J = 1, ICNT
               ISCR (IOFF + J) = 1
            END DO
            IOFF = IOFF + ICNT
         END DO

         IF (TOTSH .GT. ACESSH) THEN

            PART = MOD (TOTSH,ACESSH)
            LOOP = (TOTSH - PART) / ACESSH
            ISTART = 1
            IEND   = ISTART + ACESSH - 1
            DO I = 1, LOOP
               WRITE (UNEW,'(14I5)') (ISCR (J), J = ISTART, IEND)
               ISTART = ISTART + ACESSH
               IF ((TOTSH - IEND) .GT. ACESSH) THEN
                  IEND = IEND + ACESSH
               ELSE
                  IEND = TOTSH
               ENDIF
            END DO

            IF (PART .GT. 0) THEN
               WRITE (UNEW,'(14I5)') (ISCR (J), J = ISTART, IEND)
            END IF

         ELSE
            WRITE (UNEW,'(14I5)') (ISCR (J), J = 1, TOTSH)
         END IF
C
C
C             ...Now write out the the number of primitive functions
C                for that contracted orbital of given angular momenta!
C
C
         IOFF = 0
         DO I = 1, NSHELL
            ICNTR = NCNTSH (I)
            IPRIM = NPRMSH (I)
            DO J = 1, ICNTR

               ICOUNT = 0
               KK     = 0
               DO K = 1, IPRIM
                  DCHECK = DABS (COEFS (I,J,K))
                  IF (DCHECK .GT. CTHRESH) THEN
                     ICOUNT = ICOUNT + 1
                     KK     = KK + 1
                     MEXP (I,J,KK) = EXPNT (I,K)
                     MCEF (I,J,KK) = COEFS (I,J,K)
                  END IF
               END DO

               ISCR (IOFF + J) = ICOUNT

            END DO
            IOFF = IOFF + ICNTR
         END DO

         IF (TOTSH .GT. ACESSH) THEN

            PART = MOD (TOTSH,ACESSH)
            LOOP = (TOTSH - PART) / ACESSH
            ISTART = 1
            IEND   = ISTART + ACESSH - 1
            DO I = 1, LOOP
               WRITE (UNEW,'(14I5)') (ISCR (J), J = ISTART, IEND)
               ISTART = ISTART + ACESSH
               IF ((TOTSH - IEND) .GT. ACESSH) THEN
                  IEND = IEND + ACESSH
               ELSE
                  IEND = TOTSH
               ENDIF
            END DO

            IF (PART .GT. 0) THEN
               WRITE (UNEW,'(14I5)') (ISCR (J), J = ISTART, IEND)
            END IF

         ELSE
            WRITE (UNEW,'(14I5)') (ISCR (J), J = 1, TOTSH)
         END IF

         WRITE (UNEW,  *  )  ''
C
C
C             ...To summarize so far:
C
C                  ISCR contains the number of primitives per contracted
C                  function.
C
C                  MEXP contains the exponents of the contracted function
C                  all moved to the top.
C
C                  MCEF contains the coefficients matched with the exponents
C                  in MEXP.
C
C             ...All we have to do now is loop over and write the file!
C
C
         IOFF = 0
         DO I = 1, NSHELL
            ICNTR = NCNTSH (I)
            DO J = 1, ICNTR
               IPRIM = ISCR (IOFF + J)
               WRITE (UNEW,'(5F14.6)')  (MEXP (I,J,K), K = 1, IPRIM)
               WRITE (UNEW,*) ''
               DO K = 1, IPRIM
                  WRITE (UNEW,'(F10.7,1X)') MCEF (I,J,K)
               END DO
               WRITE (UNEW,*) ''
            END DO
            IOFF = IOFF + ICNTR
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END

