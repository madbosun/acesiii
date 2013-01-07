#include "flags.h"

      Subroutine A3_ECP_MAIN(Iecp, Natoms, Namat, Iqmstr, Jcostr,
     &                       Nucstr, Nrcstr, Jstrt, Jrs, Expa, 
     &                       Nalpha, Cont, Npcoef, Cord, Charge,
     &                       IGenby, Nharm)
C      
      Implicit Double Precision(A-H, O-Z)

#include "ecp.par"
#include "ecp.com"
#include "flags.com"
#include "flags2.com"

      Character*2 mtype(Max_centers)
      Character*80 Cnick(Max_centers)
      Character*80 Ecpnam(Max_centers)
      Character*4 namat(max_centers)
      Logical nharm

      Dimension iqmstr(Max_centers),jcostr(Max_centers,maxjco),
     &          nucstr(Max_centers,Maxang,Maxjco),
     &          nrcstr(Max_centers,Maxang,Maxjco),
     &          Jstrt(Max_shells),Jrs(Max_shells),Expa(Nalpha),
     &          Cont(Npcoef),Cord(3,Natoms),
     &          IGenby(Max_centers), Charge(Max_centers),
     &          Fcord(3,Max_centers)

      Dimension Ecpint(Max_cbf*Max_cbf), Ecpint_C(Max_cbf*Max_cbf)

C ECP integrals. Except for simplifications, I have kept the original
C form of all the input/output ECP processing routines. This was 
C originaly written by Christien Huber, Ajith Perera, 12/2001.
C
C - Natoms is the number of symmetry unique atoms
C - Ipseux is 0 when no ECP definitions; > 0 otherwise.
C - Cnick is the name of then ECP 
C - iqmstr highest ang. mom. qun. num.(AMQN) for each atom.
C - jcostr has the AQMN for each shell (each atom). Normally the # of shells
C   is the same as the Max AMQN but some cases there may be mutiple shells
C   with the same AQMN, so jcostr tells you how many.
C - nucstr is the number of primitives per each shell (per atom,per AQMN).
C - nrcstr is the number of contracted functions per each shell
C   (per atom,per AQMN)
C - nhram is logical that tells you whether this is sperical or 
C   cartesian basis calculation (key-word driven)
C - I must add that the names choosen for those 4 variables are dumb!!
C   
      Nharm = .False. 
      Write(6,*) "Spherical flags", Iflags(h_IFLAGS_spherical)
      If (Iflags(h_IFLAGS_spherical) .EQ. 1) Nharm = .True. 
CSS      Itol  = Iflags(h_IFLAGS2_intgrl_tol)
CSS      Tol   = Dble(10.0D+00*(-Itol))
 
      if (iecp.ne.0) then
c-----------------------------------------------------------------------
c get atomsymbols and put them to mtype-array
c-----------------------------------------------------------------------
         ierr=0
         call chrges(CHARGE,mtype,natoms,ierr)
c-----------------------------------------------------------------------
c read ecp data and prepare data for evulation of the c integrals
c-----------------------------------------------------------------------

         if (ierr.eq.0) then
             ipr=0
             call ecpdef(natoms,mtype,Ecpnam,cnick)

             if (ierr.eq.0) then
                ipr=0
                call ecppar(ipr,cnick,CHARGE,natoms,.FALSE.,0)
C
c-----------------------------------------------------------------------
c  Summation over the charges and storage of this value
c-----------------------------------------------------------------------

                chgsum=0
                do i=1,natoms
                    chgsum=chgsum+charge(i)
                Enddo

                nproton=idint(chgsum)
c
                call putrec(20,'JOBARC','NMPROTON',1,nproton)

c-----------------------------------------------------------------------
c prepare data for ecp integral evaluation
c-----------------------------------------------------------------------
                Call ecp_init
                Call A3_ecp_basis_init(Natoms, Namat, Iqmstr, Jcostr, 
     &                                 Nucstr, Nrcstr, Jstrt, Jrs, 
     &                                 Expa, Nalpha, Cont, Npcoef,
     &                                 Nharm)
c-----------------------------------------------------------------------
             else
                write(6,"(a)") 'ERROR while reading ECP data!'
                write(6,"(a)") 'Please check yourECPDATA file.'
                Call Errex
             endif

         else
             write(6,"(a)") 'E R R O R in subroutine chrges!'
             Call Errex
         endif
c-----------------------------------------------------------------------

      Endif

      call prepareecp(Cord, Charge, Igenby, Natoms, Ntotatoms)
      call A3_ecp_int_driver(Ecpint, Ecpint_C, Natoms, Ntotatoms,
     &                       Igenby, Cord, .FALSE., Nharm)

#ifdef _OBSOLETED
c-----------------------------------------------------------------------
c Test of numerical ECP derivative integrals. VFL 6/14/2012 
c-----------------------------------------------------------------------

     call set_gradpoints(max_centers,Cord, Natoms, Displ, eps) 

c-----------------------------------------------------------------------
c Save the original coordinates 
c-----------------------------------------------------------------------

      do ii = 1, Natoms 
      do jj = 1, 3 
         FCord(jj,ii) = Cord(jj,ii) 
      enddo ! jj 
      enddo ! ii 

      open(unit=62, file='decpint.dat', form='unformatted')
      rewind 62

      write(62) Natoms 
      do i = 1, Natoms 
      do j = 1, 3 

         do jj = 1, Max_cbf*Max_cbf  
            Decpint(i,j,jj)  = 0.0 
         enddo 

      do k = 1, 2 

         do ii = 1, Natoms 
         do jj = 1, 3 
             Cord(jj,ii) = FCord(jj,ii) 
         enddo ! jj 
         enddo ! ii 

         Cord(j,i) = Displ(i,j,k) 
         write(6,*) ' DISP X :', i, j, Cord(j,i) 

         do jj = 1, Max_cbf*Max_cbf  
            Ecpint(jj) = 0.0 
         enddo 

CSSS         call A3_ecp_nuder_driver(Ecpint,Cord, Charge, Natoms, 
CSSS     *                            nint, i, j, k)

         if (k .eq.1) then ! + displacent 
            do jj = 1, nint*(nint+1)/2   
               Decpint(i,j,jj)  = Ecpint(jj)  
            enddo 
         endif 
c        write(6,*) 'VFL',(ecpint(jj), jj =1, nint*(nint+1)/2)  

c        rconv = 0.529177246 
         if (k .eq.2) then ! - displacent 
            do jj = 1, nint*(nint+1)/2   
               Decpint(i,j,jj)  = (Decpint(i,j,jj) - Ecpint(jj))
     *                           /(Displ(i,j,1) - Displ(i,j,2)) ! eps    
            enddo 
         endif 

      enddo ! k 

         write(62) nint  
         write(62) (Decpint(i,j,jj), jj =1, nint*(nint+1)/2)  
c        write(6,*) 'VFL',(Decpint(i,j,jj), jj =1, nint*(nint+1)/2)  

      enddo ! j 
      enddo ! i 

      close(62) 
      write(6,*) ' ECP derivative file closed'  

c     do ii = 1, Natoms 
c     do jj = 1, 3 
c         Cord(jj,ii) = FCord(jj,ii) 
c     enddo ! jj 
c     enddo ! ii 
c     call A3_ecp_grdint_driver(Ecpint,Cord, Charge, Natoms)
#endif
      
      Return
      End
