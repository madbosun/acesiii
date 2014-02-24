      subroutine rad_mesh(r_point,r_weight, 
     & radnum,bsrd,atnumber)
           
      implicit none
      
      integer i, j, k, l,radnum,atnumber
      double precision r_point(radnum),bsrd(105)
      double precision r_weight(radnum),kkk
      
      
      
c      write(6,*) 'Atmic number    ',atnumber
c      write(6,*) 'Bregg - Slater parameter ',bsrd(atnumber)   
      
      do i=0,radnum-1
      kkk = dble(i)/dble(radnum)
      r_point(i+1) = bsrd(atnumber)*( kkk**2.d0 )/
     & ( (1.d0-kkk)**2.d0)
      r_weight(i+1) = bsrd(atnumber)*2.d0*kkk/
     & ( (1.d0-kkk)**3.d0)  
      
      enddo
        
      return
      end
c------------------------------------------------------------------     
c                New version of radial mesh                       |
c------------------------------------------------------------------
 
      subroutine rad_mesh_new(r_point,r_weight,radnum,bsrd,atnumber)
      implicit none

      integer I,radnum,atnumber
      double precision r_point(radnum),r_weight(radnum),QQQ,PPP,
     & bsrd(105)
         PPP = dble(radnum+1)
         DO I=1,radnum
         QQQ =  dble(radnum+1-I)
         r_point(I) = bsrd(atnumber)*( (dble(I))**2.d0 )/(QQQ**2.d0)
         r_weight(I) = bsrd(atnumber)**3.d0* 
     &            ( 2.d0*PPP*(dble(I))**5.d0 )/(QQQ**7.d0)
         ENDDO 

      return
      end
