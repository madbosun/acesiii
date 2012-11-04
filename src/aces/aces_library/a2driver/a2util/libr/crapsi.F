
c This routine was the ACES initialization routine, but that function has been
c transferred to aces_init. However, this routine is not merely kept for
c compatibility. It is still needed to initialize some legacy common blocks;
c thus, removing them from the acescore namespace.

c NO NEW MEMBER EXECUTABLES SHOULD CALL THIS!

c To use the automatic memory allocation, IENTRY must be >= 0.

c#define _DEBUG_CRAPSI

      SUBROUTINE CRAPSI(ICORE,IUHF,IENTRY)
      IMPLICIT INTEGER (A-Z)

      INTEGER ICORE(*), IUHF, IENTRY

      LOGICAL BTMP

c COMMON BLOCKS
      COMMON /ISTART/ I0,ICRSIZ
c   o from the old MCHPRM routine
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /BUFLEN/ ILNBUF
      LOGICAL MINTPRC,MVDINT1,MVDINT2,MANTI
      COMMON /MINTPC/ MINTPRC
      COMMON /MVDINT/ MVDINT1,MVDINT2
      COMMON /MEMANTI/ MANTI
c   o from the old GETSTF routine
      LOGICAL LINCC,CICALC
      COMMON /LINEAR/ LINCC,CICALC

c ----------------------------------------------------------------------

c   o initialize the ACES environment
      bTmp = (ientry.gt.-1)
      call aces_init(icore,i0,icrsiz,iuhf,bTmp)

c ----------------------------------------------------------------------

c   o initialize legacy data

c     /BUFLEN/ (only used in intprc)
      ILNBUF = 1536

c     /FILES/
      LUOUT  = 6
      MOINTS = 50

c     /LINEAR/
      LINCC  = .FALSE.
      CICALC = .FALSE.

#if defined(_UNICOS) || ( defined(_AIX) && defined(F_64BIT) )
      MINTPRC = .FALSE.
      MVDINT2 = .FALSE.
      MANTI   = .FALSE.
#else
      MINTPRC = .TRUE.
      MVDINT2 = .TRUE.
      MANTI   = .TRUE.
#endif /* _UNICOS || (_AIX && F_64BIT) */

c ----------------------------------------------------------------------

      RETURN
      END

