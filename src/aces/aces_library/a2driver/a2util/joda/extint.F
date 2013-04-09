      Subroutine ExtInt (XTol, NrX, AtNrX, QX, NrI, AtNrI, QI,
     $   XtoI, IStat)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose: Determine the set of pointers (XtoI) which will take a center
C          in the external representation of the list of centers and
C          "point to" the correct center in the internal list.
C
C Modified 7/91 to write pointer list to JOBARC by JFS
C
C Arguments:
C     NrX     Number of external centers (input)
C     AtNrX   Atomic number list for the external set (input)
C     QX      Cartesian coordinates for the external set (input)
C     NrI     Number of internal centers (input)
C     AtNrI   Atomic number list for the internal set (input)
C     QI      Cartesian coordinates for the internal set (input)
C     Tol     Accuracy required of "identical" coordinates (input)
C     XtoI    External to internal permutation list (output)
C             Each element gives the internal rep. center number
C             corresponding to the external rep. center number.
C     IStat    Returns success or error in finding XtoI (output)
C             = 0  Successfully computed permutation list.
C             = 1  Mismatch in number of (non-dummy) centers.
C             = 3  Mismatch in coordinates.
C Scratch arguments:
C     Scr     Space to work on coordinates (scratch)
C             (Room for 3 sets of cartesian coordinates)
C     OrdX    Work space for permutation list (scratch)
C     OrdI    Work space for permutation list (scratch)
C
C Dependents:
C     SortC   Sort cartesian coordinates, weight with at. number
C     SAXPY   SCIPORT BLAS routine (double precision, despite the name)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C $Log: extint.F,v $
C Revision 1.1.1.1  2003/04/02 19:21:35  aces
C INITIAL 2.4 IMPORT
C
C Revision 4.0  89/03/14  01:14:18  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C
C Revision 3.5  89/02/16  14:32:44  bernhold
C Initial version
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer NrX, AtNrX(NrX)
      Integer NrI, AtNrI(NrI)
      Integer XtoI(NrX), IStat
      Double precision QX(3*NrX), QI(3*NrI)
      Double precision Tol,XTol
C
C     Scratch arguments - should eventually use ED's dynamic alloc.
C
#include "mxatms.par"
      PARAMETER (TOL=1.D-4)
      Integer OrdX(MxAtms), OrdI(MxAtms)
      Double precision Scr(3*3*MxAtms)
C
C     Define the "Standard I/O" units to use throughout the code
      Parameter (LuIn  = 5)
      Parameter (LuOut = 6)
      Parameter (LuErr = 6)
C
      Integer i
C
C     Sort both sets of coordinates
C
      Call SortC (NrX, QX, AtNrX, Scr(1),       OrdX, Scr(6*NrI+1) )
      Call SortC (NrI, QI, AtNrI, Scr(3*NrI+1), OrdI, Scr(6*NrI+1) )
C
C     Check to make sure they're the same coordinates!
C     First, be sure that the only things missing from X are dummies.
C
      Do 10 i = NrX+1, NrI
         If ( AtNrI(OrdI(i)) .ne. 0) then
            IStat = 1
            Write (LuErr, 9910)
            Return
         EndIf
 10   Continue
 9910 Format (' @EXTINT-F, Number of internal and external centers ',
     $   'do not match.')
C
C     Second, be sure the coordinates themselves match.
C
      Call SAXPY (3*NrX, -1.0d0, Scr(3*NrI+1), 1, Scr(1), 1)
      Do 20 i = 1, 3*NrX
         DIFF=SCR(I)
         If (DIFF.GT.TOL)THEN
            IStat = 3
            Write (LuErr, 9940)
            Write (LuErr, 9935) (j, AtNrX(j), QX(3*(j-1)+1),
     $         QX(3*(j-1)+2), QX(3*(j-1)+3), j = 1, NrX)
            Write (LuErr, 9945)
            Write (LuErr, 9935) (j, AtNrI(j), QI(3*(j-1)+1),
     $         QI(3*(j-1)+2), QI(3*(j-1)+3), j = 1, NrI)
            Write (LuErr, 9930)
            Return
         EndIf
 20   Continue
 9930 Format (' @EXTINT-F, External and internal cartesian',
     $   ' coordinates do not match.')
 9935 Format (' Center',3X,'At. Nr.',3X,19X,'X',19X,'Y',19X,'Z'/
     $   (I7,3X,I7,3X,3F20.10))
 9940 Format (/' External cartesian coordinates:')
 9945 Format (/' Internal cartesian coordinates:')
C
C     Now make the permutation list
C
      Do 100 i = 1, NrX
         XtoI(OrdX(i)) = OrdI(i)
 100  Continue
C
      CALL PUTREC(20,'JOBARC','EXTINTVC',NRX,XTOI)
      Return
      End
