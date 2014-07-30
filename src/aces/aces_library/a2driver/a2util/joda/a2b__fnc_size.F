         SUBROUTINE  A2B__FNC_SIZE
     +
     +                    ( UACES,NSHELL,
     +                      MAXCNTR,MAXPRIM,
     +
     +                                 LANGSH,
     +                                 NCNTSH,NPRMSH )
     +
     +
C-----------------------------------------------------------------------
C  OPERATION   : A2B__FNC_SIZE
C  MODULE      : ACESII BASIS CONVERTER
C  MODULE-ID   : A2B
C  SUBROUTINES : none
C  DESCRIPTION : This routine will determine the number of contracted
C                functions and primitive functions that there are.
C
C                Input:
C
C                     UACES   = The unit number for the ZMAT.BAS
C                     NSHELL  = The number of shells in the ZMAT.BAS
C                     MAXCNTR = The maximum number of contracted
C                               functions in the basis set
C                     MAXPRIM = The maximum number of primitive 
C                               functions in the basis set
C
C                Output:
C
C                     LANGSH  = An array of angular momentum numbers
C                     NCNTSH  = An array of the number of contracted
C                               functions for a given angular momenta
C                     NPRMSH  = An array of the number of primitive
C                               functions for a given angular momenta
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

         INTEGER    JUNK
         INTEGER    I
         INTEGER    UACES
         INTEGER    NSHELL,MAXCNTR,MAXPRIM

         INTEGER    LANGSH (1:NSHELL)
         INTEGER    NCNTSH (1:NSHELL)
         INTEGER    NPRMSH (1:NSHELL)
C
C
C------------------------------------------------------------------------
C
C
C             ...Read the size of functions in the ZMAT.BAS
C
C        
         READ (UACES,'(14I5)') (LANGSH (I), I = 1, NSHELL)
         READ (UACES,'(14I5)') (NCNTSH (I), I = 1, NSHELL)
         READ (UACES,'(14I5)') (NPRMSH (I), I = 1, NSHELL)
         READ (UACES,'(A)') JUNK

         MAXCNTR = MAXVAL (NCNTSH)
         MAXPRIM = MAXVAL (NPRMSH)
C
C
C             ...ready!
C
C
         RETURN
         END

