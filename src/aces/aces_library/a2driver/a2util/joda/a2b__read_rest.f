         SUBROUTINE  A2B__READ_REST
     +
     +                    ( UACES,
     +                      NSHELL,MAXCNTR,MAXPRIM,
     +                      LANGSH,NCNTSH,NPRMSH,
     +                      SCR,
     +
     +                                     EXPNT, COEFS )
     +
     +
C-----------------------------------------------------------------------
C  OPERATION   : A2B__READ_REST
C  MODULE      : ACESII BASIS CONVERTER
C  MODULE-ID   : A2B
C  SUBROUTINES : none
C  DESCRIPTION : This routine will read the exponents and coefficients
C                of a basis set.
C
C                Input:
C
C                     UACES   = The unit number for the ZMAT.BAS
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
C                     SCR     = A double precision scratch array
C
C                Output:
C
C                     EXPNT   = A matrix of the exponents for a given
C                               angular momenta
C                     COEFS   = A matrix of the coefficients for a 
C                               given angular momenta and contracted
C                               function
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

         INTEGER    I,J,K,L,JUNK
         INTEGER    ISTART,IEND,PART,LOOP
         INTEGER    IPRIM,ICNTR
         INTEGER    UACES,ACESEXP,ACESCNT
         INTEGER    NSHELL,MAXCNTR,MAXPRIM

         DOUBLE PRECISION  ZERO

         INTEGER    LANGSH (1:NSHELL)
         INTEGER    NCNTSH (1:NSHELL)
         INTEGER    NPRMSH (1:NSHELL)

         DOUBLE PRECISION  EXPNT (1:NSHELL,1:MAXPRIM)
         DOUBLE PRECISION  COEFS (1:NSHELL,1:MAXCNTR,1:MAXPRIM)
         DOUBLE PRECISION  SCR   (1:NSHELL,1:MAXPRIM,1:MAXCNTR)

         PARAMETER  ( ZERO    = 0.0D0 )
         PARAMETER  ( ACESEXP = 5     )
         PARAMETER  ( ACESCNT = 7     )
C
C
C------------------------------------------------------------------------
C
C
C             ...Zero out the arrays first!
C
C
         DO I = 1, NSHELL
         DO J = 1, MAXPRIM

            EXPNT (I,J) = ZERO

            DO K = 1, MAXCNTR

               SCR   (I,J,K) = ZERO
               COEFS (I,K,J) = ZERO

            END DO
         END DO
         END DO
C
C
C             ...Read the remainder of the ZMAT.BAS
C
C
         DO I = 1, NSHELL

            IPRIM = NPRMSH (I)
            ICNTR = NCNTSH (I)
C
C
C             ...Read in the exponents for the given shell!
C
C
            IF (IPRIM .GT. ACESEXP) THEN

               PART = MOD (IPRIM,ACESEXP)
               LOOP = (IPRIM - PART) / ACESEXP
               ISTART = 1
               IEND   = ISTART + ACESEXP - 1
               DO J = 1, LOOP
                  READ(UACES,'(5F14.7)') (EXPNT (I,K), K = ISTART, IEND)
                  ISTART = ISTART + ACESEXP
                  IF ((IPRIM - IEND) .GT. ACESEXP) THEN
                     IEND = IEND + ACESEXP
                  ELSE
                     IEND = IPRIM
                  END IF
               END DO
            
               IF (PART .GT. 0) THEN
                  READ(UACES,'(5F14.7)') (EXPNT (I,K), K = ISTART, IEND)
               END IF

            ELSE
               READ(UACES,'(5F14.7)') (EXPNT (I,K), K = 1, IPRIM)
            END IF

            READ (UACES,'(A)') JUNK
C
C
C             ...Read in the contraction coefficients for the given
C                shell and primitive!
C
C
            DO J = 1, IPRIM

               IF (ICNTR .GT. ACESCNT) THEN

                  PART = MOD (ICNTR,ACESCNT)
                  LOOP = (ICNTR - PART) / ACESCNT
                  ISTART = 1
                  IEND   = ISTART + ACESCNT - 1
                  DO L = 1, LOOP

                     READ (UACES,'(7(F10.7,1X))')
     +
     +                    (SCR (I,J,K), K = ISTART, IEND)
     +
                     ISTART = ISTART + ACESCNT
                     IF ((ICNTR - IEND) .GT. ACESCNT) THEN
                        IEND = IEND + ACESCNT
                     ELSE
                        IEND = ICNTR
                     ENDIF
                  END DO

                  IF (PART .GT. 0) THEN

                     READ (UACES,'(7(F10.7,1X))')
     +
     +                    (SCR (I,J,K), K = ISTART, IEND)
     +
                  END IF

               ELSE

                     READ (UACES,'(7(F10.7,1X))')
     +
     +                    (SCR (I,J,K), K = 1, ICNTR)
     +
               END IF

            END DO

            READ (UACES,'(A)') JUNK
C
C
C             ...End the loop over shells!
C
C
         END DO
C
C
C             ...We need to transpose the SCR array so that the
C                coefficient array lines up in a better order!
C
C
         DO I = 1, NSHELL
         DO J = 1, MAXCNTR
         DO K = 1, MAXPRIM

            COEFS (I,J,K) = SCR (I,K,J)

         END DO
         END DO
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END

