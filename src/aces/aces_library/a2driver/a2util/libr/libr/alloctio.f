      Subroutine AllocTIO(MAvail, Reserve)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: alloctio.f,v 1.1.1.1 2003/04/02 19:21:39 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     AllocTIO -- Determine memory to reserve for temporary I/O vectors
C
C SYNOPSIS
C&      Include 'implicitnone.inc'
      Integer MAvail, Reserve
C
C ARGUMENTS
C     MAvail  Words of memory available after all permanent I/O needs
C             are satisfied [IN]
C     Reserve Number of words of memory required for temporary I/O
C             vectors [OUT]
C
C DESCRIPTION
C     Estimate the amount of memory required for any temporary I/O
C     vectors we might generate later in the code.  Decision may be
C     based on the particular model being calculated, the memory
C     available, etc.
C
C COMMON BLOCKS
C&      Include 't3ctrl.inc'
C
C     Control execution of triples contributions
C     T3Meth(1)  Level of T3 inclusion.
C                < 0  Non-iterative models:
C                     -3 => FS CCSD+T(3)
C                     -4 => FS CCSD+T(CCSD) [not yet implemented]
C                = 0  No triples.
C                > 0  Iterative triples models [not yet implemented]
C     T3Meth(2)  Symmetries of T3(0,1)(lab,kij) required.
C                = 1  Only totally symmetric (k,l).
C                = 8  All symmetries of (k,l). [not yet implemented]
C     T3Meth(3)  Types of T3(0,1)(lab,kij) required.
C                = 0  Spectator amplitudes (k=l).
C                = 1  Active-active block of (k,l).
C                = 2  All possible (k,l).
C     T3Meth(4)  Use H or Hbar in contractions
C                = 0  Use H
C                = 1  Use Hbar
C
C     See T3MetInf for information about the "standard" (named) triples
C     models in their T3Meth vectors.
C
      Integer LT3Meth
      Parameter (LT3Meth = 4)
      Integer T3Meth(LT3Meth)
      Common /T3Ctrl/ T3Meth
C
      Integer NOccO(2), NVrtO(2)
      Common /Info/ NOccO, NVrtO
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     I       General Use.
C     J       General Use.
C
      Integer I, J
C
C     Basic CCSD requires no temporary space.
C
      Reserve = 0
C
C     Handle T3 models which compute vectors for each set of fixed
C     indices: spectator-T(3), T(3)
C
C     There are 2 fixed indices (in general, 1 for spectator-only models).
C     There are 11 relevant ITypes for each spin of the fixed indices.
C     The longest is at most MAX( virtual, occupied) orbitals.
C     Four of the ITypes are UPLT or UPLE, which require a "zero list"
C     as long as their gather vector.
C
C     *** If accomodating multiple (K,L) in a single pass ***
C     The number of fixed indices done at a time is bounded by the number
C     of occupied orbitals in one case and number of active in the other.
C     Allocate for the number of occupied for both -- very generous.
C
      I = Max( NOccO(1), NOccO(2) )
      J = Max( I, NVrtO(1), NVrtO(2) )
C
      If (T3Meth(1) .eq. -3) then
         If (T3Meth(3) .eq. 0) then
            Reserve = (11 + 6) * J
         Else
            Reserve = 2 * (11 + 6) * J
         EndIf
      EndIf
C
      Return
      End
