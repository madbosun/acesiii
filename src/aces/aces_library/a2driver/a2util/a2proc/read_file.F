       Subroutine Read_file(NATOMS, NVIBS, Input_Array, Scratch1)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "flags.h"
#include "machsp.com"
#include "jodaflags.com"                   
#include "mxatms.par"
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
      

#ifdef _DEBUG_LVL0
      Write(6,*) "Nfull and Nleft", Nfull, Nleft
#endif

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

#ifdef _DEBUG_LVL0
      nt = nfull + Nleft
      Write(6,*) "Data as read from the external file"
      Call output(Scratch1, 1, Natoms*nt, 1, 9, natoms*(nfull+3), 9, 
     &           1)
#endif
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

#ifdef _DEBUG_LVL0
      nt = nfull + Nleft
      Write(6,*) "Reordered data from  external file"
      Call output(Input_array, 1, Natoms*3, 1, natoms*(nfull+3), 9, 
     &           1)
#endif
      Return
      End
