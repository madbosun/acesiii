
C FILLS OUT MACHINE-DEPENDENT COMMON BLOCKS

      SUBROUTINE MCHPRM

      LOGICAL MVSCF,MINTPRC,MVDINT1,MVDINT2,MANTI
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /MACHSP2/MASK1,MASK2,ISHFSZ
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FILSPC/ ILNBUF,IPRCLN,IPRCWD
      COMMON /IOPOS/ ICRSIZ,ICHCSZ,IOFF(2),LENREC
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MEMSCF/ MVSCF
      COMMON /MINTPC/ MINTPRC
      COMMON /MVDINT/ MVDINT1,MVDINT2
      COMMON /MEMANTI/ MANTI
#include "icdacc.com"

C IINTLN IS THE INTEGER LENGTH IN BYTES.
C IFLTLN IS THE FLOATING POINT LENGTH IN BYTES.
C IINTFP IS THE RATIO IFLTLN/IINTLN

#ifdef F_64BIT
      IINTLN=8
#else
      IINTLN=4
#endif
      IFLTLN=8
      IINTFP=(IFLTLN/IINTLN)
      IALONE=2**(8*IINTLN/4)-1
      IBITWD=8*IINTLN/4

cYAU - Okay, what is this "idacc?" all about...
c
c    First off, the name refers to Direct ACCess. The ultimate problem
c is that the Fortran standard lets RECL take any unit when specifying
c the record length of direct access, unformatted files. It's not
c necessarily compiler, OS, or architecture dependent. In fact, some
c configurations let you specify the units of RECL with an environment
c variable!
c    To aid portability, idaccm was introduced by Nevin as a scaling
c factor for RECL. In essence, idaccm is the number of RECL units
c per integer; therefore, it could take any of the values: 1, 2, 4, or 8.
c    If there is an accompanying Fortran 9x compiler, then a short
c program will tell you whether to define _RECL_IS_WORDS_:
c       program test
c          integer :: i, len
c          inquire(iolength=len) i
c          print *, len
c       end program test
c Whatever number is printed, that should be the value of idaccm. For
c the time being, if RECL is not in bytes, then we will assume 4-byte
c words. It is possible RECL is 8-byte words, but that is very rare.

#ifdef _RECL_IS_WORDS_
#ifdef F_64BIT
      idaccm=2
#else
      idaccm=1
#endif /* F_64BIT */
#else
      idaccm=iintln
#endif /* _RECL_IS_WORDS_ */

#ifdef _DEBUG
c     We will try to print this info just once. Since ACES needs files
c     to get things done, we will check for a characteristic file and
c     print only if it does not exist. Also, since MVSCF gets set a
c     few lines down, we will use it temporarily.
      inquire(file='JOBARC',exist=MVSCF)
      if (.not.MVSCF) then
         write(*,*) '@MCHPRM: Machine/Build-specific parameters:'
         write(*,*) '         bytes per integer: ', iintln
         write(*,*) '         bytes per double:  ', ifltln
         write(*,*) '         recls per integer: ', idaccm
      end if
#endif /* _DEBUG */

      ILNBUF=1800
      LUOUT=6
      MOINTS=50

      IJUNK=4*IBITWD-3
      ISHFSZ=IJUNK
      MASK1=2**IJUNK-1
      MASK2=7

C LOGICAL VARIABLES TO CONTROL MEMORY IN SEVERAL MODULES
C NOTE, ON A RISC IT IS ADVANTAGEOUS TO REDUCE MEMORY
C REQUIREMENTS IN ALL SORT STEPS.

#if (defined _UNICOS) || ( (defined _AIX) && (defined F_64BIT) )
      MVSCF=.FALSE.
      MINTPRC=.FALSE.
      MVDINT1=.TRUE.
      MVDINT2=.FALSE.
      MANTI=.FALSE.
#else
      MVSCF=.FALSE.
      MINTPRC=.TRUE.
      MVDINT1=.TRUE.
      MVDINT2=.TRUE.
      MANTI=.TRUE.
#endif /* _UNICOS || (_AIX && F_64BIT) */

      RETURN
      END

