       Subroutine Read_file(NATOMS, NVIBS, Input_Array, Scratch1)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C












































































































































































































































































































































































































































































































































c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end










c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C#include "istart.com"                      
C    
      Double Precision Input_Array
      Character*80 junk
      Character*2 Zsym(Natoms)
      Logical NORMCO_EXSIST
      Dimension Input_array(3*Nvibs,3*Nvibs), Scratch1
     &          (Natoms*(Nvibs/3+3),9)
C
      DATA XTANG /0.529177249D0/ 
C
      Iunit = 4
      INQUIRE(FILE='NORMCO',EXIST=NORMCO_EXSIST)
      IF(NORMCO_EXSIST)THEN
        OPEN(UNIT=Iunit,FILE='NORMCO',FORM='FORMATTED',STATUS='OLD')
        REWIND(Iunit)
      ELSE
        WRITE(6,20)
 20      format(T3,'@-READ_file, NORMCO file not present.')
        Call Errex
      ENDIF
     
      Ncols = 3
      Nfull = Nvibs/3
      Nleft = Nvibs - 3*Nfull
      

      Write(6,*) "Nfull and Nleft", Nfull, Nleft

      Ioff = 1
      Do Jcount = 1, Nfull
     
         Do Ijunk = 1, 4
            Read(Iunit,*) junk
         Enddo
         Do Iatoms = 1, Natoms
            Read(Iunit,10) Zsym(Iatoms), (Scratch1(Ioff, joff), 
     &                     joff=1,9)
            Ioff = Ioff + 1
         Enddo
      Enddo
  
      Do Jcount = 1, Nleft

         Do Ijunk = 1, 4
            Read(Iunit,*) junk
         Enddo
         Do Iatoms = 1, Natoms
            Read(Iunit,10) Zsym(Iatoms), (Scratch1(Ioff, joff), 
     &                     joff=1,3*Nleft)
            Ioff = Ioff + 1
         Enddo
      Enddo

      nt = nfull + Nleft
      Write(6,*) "Data as read from the external file"
      Call output(Scratch1, 1, Natoms*nt, 1, 9, natoms*(nfull+3), 9, 
     &           1)
 10   FORMAT(A2, 3X, 3(3(F7.4,1X), 2x))
	
      Ntotal = Nfull + Nleft
      Do Indx = 1, Natoms*Ntotal 
         Ioff = 0

          Do Jndx = 1, NvibS, 3
             Ioff = Ioff + 1
             
             Do ixyz =Jndx, Jndx + 2
       
                Input_array(Jndx, Ioff) = Scratch1(Indx, Jndx) 

             Enddo
          Enddo
      Enddo

      nt = nfull + Nleft
      Write(6,*) "Reordered data from  external file"
      Call output(Input_array, 1, Natoms*3, 1, natoms*(nfull+3), 9, 
     &           1)
      Return
      End
