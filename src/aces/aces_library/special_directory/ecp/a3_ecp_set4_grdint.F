#include "flags.h"
C
      Subroutine A3_ECP_set4_grdint(Iecp,Nmax,namat,iqmstr,jcostr,
     &                              nucstr, nrcstr, Jstrt, Jrs, 
     &                              Expa, Nalpha, Cont, Npcoef, 
     &                              Coord, Charge)
C      
      Implicit Double Precision(A-H, O-Z)

#include "ecp.par"
#include "ecp.com"
#include "flags.com"
#include "flags2.com"
C
      Character*2 mtype(Max_centers)
      Character*80 Cnick(Max_centers)
      Character*80 Ecpnam(Max_centers)
      Character*4 namat(max_centers)
      Logical Nharm, Grads
C
      Dimension iqmstr(Max_centers),jcostr(Max_centers,maxjco),
     &          nucstr(Max_centers,Maxang,Maxjco),
     &          nrcstr(Max_centers,Maxang,Maxjco),
     &          Jstrt(Max_shells),Jrs(Max_shells),Expa(Nalpha),
     &          Cont(Npcoef),Coord(3,Max_centers),
     &          IGenby(Max_centers), Charge(Max_centers)
C
C ECP deriavtive integrals. Except for simplifications, I have kept the 
C original  form of all the input/output ECP processing routines. This was 
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
      nharm = (Iflags(h_IFLAGS_spherical) .EQ. 1)
      Itol  = Iflags(h_IFLAGS2_intgrl_tol)
      Thres = Dble(10.0D+00*(-Itol))
 
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
c prepare data for ecp integral evaluation
c-----------------------------------------------------------------------
                Call ecp_init
                Call A3_ecp_basis_init_4grd(Natoms, Namat,
     &                                      Iqmstr, Jcostr, Nucstr, 
     &                                      Nrcstr, Jstrt, Jrs, Expa,
     &                                      Nalpha, Cont, Npcoef, 
     &                                      Nharm)
          
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

      Return
      End
