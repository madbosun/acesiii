      subroutine prepareecp(Coord, Charge, Idexatom, Natoms,
     &                      Ntotatoms)

      implicit Double Precision (a-h,o-z)                                      

#include "ecp.par"
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

CSSS      call putrec(20,'JOBARC','ATOMINDX', ntotatoms, idexatom)       

      return                                                         
      end                                                            
