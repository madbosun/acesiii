      subroutine prepareecp(Coord, Charge, Idexatom, Natoms,
     &                      Ntotatoms)

      implicit Double Precision (a-h,o-z)                                      

#include "ecp.par"
#include "ecp.com"
#include "machsp.com"

      Dimension Coord(3, Max_centers), Charge(Max_centers),
     &          idexatom(Max_centers)
      
      double precision :: dist, x(3)                                 

      call getrec(20,'JOBARC','NATOMS  ', 1, ntotatoms)           
cSSS     call getrec(20,'JOBARC','COORD   ', 3*ntotatoms*iintfp, acoord)                 
      do i = 1, natoms 
         idexatom(i) = 0                                      
      enddo 
C
#ifdef _DEBUG_LVL0
      write(6,*) ' Number of atoms :', natoms 
      write(6,*) ' Total Number of atoms :', ntotatoms  
      write(6,*) ' ACOORD :', (coord(1,ja), ja=1,natoms)  
      write(6,*) ' ACOORD :', (coord(2,ja), ja=1,natoms)  
      write(6,*) ' ACOORD :', (coord(3,ja), ja=1,natoms)  
#endif

      do ia = 1, ntotatoms                                           
         do ja = 1, natoms                                           
            dist = 0.d0                                           
            do id = 1, 3                                          
               x(id) = coord(id, ja)
               dist = dist + (x(id)-coord(id,ia))**2             
            enddo                                                 

            dist = sqrt(dist)                                     
            if (dist<1.d-6) then                                  
                idexatom(ia) = ja                                  
                goto 100                                           
            endif                                                 

         enddo                                                       
 100     continue                                                    
      enddo                                                          
C--------------
C This block of code and the ishuffle common block in ecp.com 
C was created to use in the cartesian to  spherical transformation
C but it was not needed. Ajith Perera, 12/2012.
C 
#ifdef _NODISCARD

      Ioff    = 1
      Inao(1) = 1
      Do Icent = 1, Ntotatoms
         Do Ishell = 1, Nshell
            Iatom  = katom(ishell)
            If (Idexatom(Icent) .EQ. Iatom) Then
               Idegen = Ndegen(Ishell)
               Numcoi = Numcon(Ishell)
               Ioff   = Ioff + 1
               Inao(Ioff) = Inao(ioff-1) + Idegen*Numcoi
            Endif
         Enddo
      Enddo
C       
      Do Icent = 1, Ntotatoms
         Do Ishell = 1, Nshell
            Iatom  = katom(ishell)
            If (Idexatom(Icent) .EQ. Iatom) Then
               Idegen = Ndegen(Ishell)
               Numcoi = Numcon(Ishell)
               If (Numcoi .ne. 0) Then
                   Ioff = Inao(Ishell)
                   Do Ico = 1, Numcoi
CSSS                      Ifrom = Inao(Ishell) - 1 + Ico
                      Ifrom = Ico
                      Do Ixyz = 1, Idegen 
                         Iorder(Ioff) = Ifrom
                         Ifrom        = Ifrom + Numcoi
                         Ioff         = Ioff + 1
                      Enddo
                   Enddo
               Endif
            Endif
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The ang. momentum offset and reordering array"
      Write(6,"(8(1x,i3))") (INao(i), i=1, Nshell)
      Write(6,*)
      Write(6,"(8(1x,i3))") (Iorder(i), i=1, Nbf)
      Write(6,*)
#endif
#endif

CSSS      call putrec(20,'JOBARC','ATOMINDX', ntotatoms, idexatom)       

      return                                                         
      end                                                            
