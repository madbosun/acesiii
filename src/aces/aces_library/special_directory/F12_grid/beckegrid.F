      subroutine beckegrid(RADNUM,ANGNUM,r_point,r_weight,bsrd,
     & x_lebed,y_lebed,z_lebed,lebed_weight)
      implicit none
c----------------------------------------------------------------------------
c The purpose of this subroutine is to set up grid for the integration of   |
c      multicenter objects using Becke's method of fuzzy cells              |
c                                                                           |
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'dbugcom.h'
      include 'trace.h'
      include 'dropmo.h'
      include 'sial_config_params.h'
      include 'F12_gr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
c-----------------------------------------------------------------------------


c-----------Input parameters--------------------------------------------------
      integer RADNUM,atnumber,ANGNUM,N
c----------------------------------------------------------------------------- 
      
      double precision bsrd(105),r_point(RADNUM),r_weight(RADNUM),
     & x_lebed(ANGNUM),y_lebed(ANGNUM),z_lebed(ANGNUM),
     & lebed_weight(ANGNUM),dref,GRIDX,GRIDY,GRIDZ,GRIDW,AtCoordX,
     & AtCoordY,AtCoordZ,CUTLIM,ATWT

      integer coretop,i000,natm,natmx,comppopv,compmemb,atmchr,atmcrd,
     & atmmas,I,J,K,iref,compnorb,iGRDPTX,iGRDPTY,iGRDPTZ,iGRDPTWT,
     & BLSIZE,ind,funit,ind1,ind2,AIJ,VxIJ,VyIJ,VzIJ,iRIJ,iRp,iWWT,NITR

  
      integer*8 itop,itop1,itop2,itop3
      integer NEED,NEED1,ierr,iq(3)

      parameter (BLSIZE = 500)
      parameter (funit = 34 )
      parameter (NITR = 4)

C------------------------------------------------------------------------
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (iptr, ibuf)
#else
      equivalence (ibuf, dbuf)
      common dbuf    ! blank common for scratch.
#endif

#ifdef ALTIX
      dptr = dshptr     ! set dbuf pointer to shared mem.
      iptr = ishptr     ! set ibuf pointer to shared mem.
#endif


c----------------------------------------------------------------------------
c-----------Open file to store grid points and weights-----------------------

c      open(unit=funit,file='GRIDDATA',FORM='UNFORMATTED',
c     &     STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      
c-----------Reading information about atomic centers----------------------------------
      call getrec (20, 'JOBARC', 'NATOMS  ', 1, natmx)
      call getrec (20, 'JOBARC', 'NREALATM', 1, natm)
      call getrec (20, 'JOBARC', 'COMPNORB', 1, compnorb)

c------Memory allocation and partitioning---------------------------------------------
      comppopv = 1
      compmemb = comppopv + natmx
      atmchr   = compmemb + natmx
      NEED1 = atmchr + natmx
       
    
      call mem_alloc(ibuf, NEED1, intsize, itop1, .true., ierr)
      if(ierr.ne.0) then
      write(6,*) 'Not enough memory!!!!'
      return
      endif 

C-------------------------------------------------------------------------------------

      call igetrec (20, 'JOBARC', 'COMPPOPV', 
     &                            natmx, ibuf (itop1+comppopv-1))
      call igetrec (20, 'JOBARC', 'COMPMEMB', 
     &                            natmx, ibuf (itop1+compmemb-1))
      call igetrec (20, 'JOBARC', 'ATOMCHRG', 
     &                              natmx, ibuf(itop1+atmchr-1)  )


c--------------------------------------------------------------------------------------
      atmcrd  = 1
      atmmas   = atmcrd + natmx * 3
      NEED = atmmas + natmx

      call mem_alloc(dbuf, NEED, bytes_per_double, itop, .true., ierr)
      if(ierr.ne.0) then
      write(6,*) 'Not enough memory!!!!'
      return
      endif


      call getrec (20, 'JOBARC','COORD   ',natmx*3,
     &                                    dbuf(itop+atmcrd-1) )
      call getrec (20, 'JOBARC','ATOMMASS',natmx,  
     &                                    dbuf(itop+atmmas-1) )
c

      call rmdummy (
     &     natmx, natm, compnorb, ibuf (itop1+comppopv-1),
     &     ibuf(itop1+compmemb-1),ibuf (itop1+atmchr-1), 
     &     dbuf(itop+atmmas-1), dbuf(itop+atmcrd-1) )
c---------------------------------------------------------------------------------------


        call recharge(natmx, ibuf (itop1+atmchr-1) )

           
c--Step 1 - set up Slater- Bragg radii for atoms and construction of AIJ paramters-------
c----and interatomic pair vectors for the construction of cutoff functions---------------

       AIJ = 1
       VxIJ = AIJ + natmx*natmx
       VyIJ = VxIJ + natmx*natmx
       VzIJ = VyIJ + natmx*natmx
       iRIJ = VzIJ + natmx*natmx
       iRp = iRIJ + natmx*natmx
       iWWT = iRp + natmx 
       NEED = iWWT + natmx

      call mem_alloc(dbuf, NEED, bytes_per_double, itop2, .true., ierr)
      if(ierr.ne.0) then
      write(6,*) 'Not enough memory!!!!'
      return
      endif


      call bsr(bsrd)
      
      call make_aij(dbuf(itop2+AIJ-1),dbuf(itop2+VxIJ-1),
     & dbuf(itop2+VyIJ-1),dbuf(itop2+VzIJ-1),dbuf(itop2+iRIJ-1),
     & ibuf (itop1+atmchr-1),dbuf(itop+atmcrd-1),bsrd,natmx) 
       
       
c-------Step 2 - set up angular Lebedev grid-------------------------
 
      if(ANGNUM.eq.194) then
       call LD0194( x_lebed,y_lebed,z_lebed,lebed_weight,N )
      endif 


c---------------------------------------------------------------------

       ind = 0
       ind1 = 0  
       CUTLIM = 1.d-14 
                  DO K=1,natmx          

                  atnumber = iref (K, ibuf (itop1+atmchr-1))
                  AtCoordX = dref (1+3*(K-1), dbuf(itop+atmcrd-1) )
                  AtCoordY = dref (2+3*(K-1), dbuf(itop+atmcrd-1) )
                  AtCoordZ = dref (3+3*(K-1), dbuf(itop+atmcrd-1) )

c-------Step 3 - set up radilal mesh for each atom--------------------
       call zero(r_point, RADNUM)
       call zero(r_weight,RADNUM)       

       call rad_mesh_new(r_point,r_weight,RADNUM,bsrd,atnumber)

c-------Construction of the grid for each atom------------------------

          DO I=1,RADNUM
           DO J=1,N
           

           GRIDX = r_point(I)*x_lebed(J)      
           GRIDY = r_point(I)*y_lebed(J) 
           GRIDZ = r_point(I)*z_lebed(J) 

c---Here we create weights for each poit using Becke fuzzy cell method-


c---------------------------------------------------------------------- 
         call atomweights(K,natmx,dbuf(itop2+AIJ-1),dbuf(itop2+VxIJ-1),
     & dbuf(itop2+VyIJ-1),dbuf(itop2+VzIJ-1),dbuf(itop2+iRIJ-1),GRIDX,
     & GRIDY,GRIDZ,dbuf(itop2+iRp-1),dbuf(itop2+iWWT-1),ATWT,NITR)

           GRIDX = r_point(I)*x_lebed(J) + AtCoordX
           GRIDY = r_point(I)*y_lebed(J) + AtCoordY
           GRIDZ = r_point(I)*z_lebed(J) + AtCoordZ


           GRIDW = (4.d0*3.1415926535898d0)*
     &        r_weight(I)*lebed_weight(J)

c           write(6,*) r_weight(I),I,'<<----'

           GRIDW = GRIDW*ATWT


c---------------Processing of the grid---------------------------------
           call procgrid(GRIDX,GRIDY,GRIDZ,GRIDW,BLSIZE,ind,ind1,ind2,
     &     funit,CUTLIM)
c---------------------------------------------------------------------
           
           ENDDO
          ENDDO


                  ENDDO

c------------Write on disk remaining point and weights ---------------
c           if(ind2.eq.1) then
c           CALL Gwrite(funit,dbuf(itop3+iGRDPTX-1),BLSIZE  )
c           CALL Gwrite(funit,dbuf(itop3+iGRDPTY-1),BLSIZE  )
c           CALL Gwrite(funit,dbuf(itop3+iGRDPTZ-1),BLSIZE  )
c           CALL Gwrite(funit,dbuf(itop3+iGRDPTWT-1),BLSIZE )
c           endif


        write(6,*) 'Number of points on Lebedev grid ... ',N
        write(6,*) 'Number of points on radial grid ... ',RADNUM
        write(6,*) 'Number of grid points per atom ...',N*RADNUM
        write(6,*) 'Total number of grid points ...',ind1
        NPT = ind1
c      close(funit,status='keep')       
 
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine atomweights(K,natmx,AIJ,VxIJ,VyIJ,VzIJ,RIJ,GRIDX,GRIDY,
     & GRIDZ,Rp,WWT,ATWT,NITR)
      implicit none
c-----------------------------------------------------------------------
c This subroutine compute atomic weights for each grid point           |
c  using Becke fuzzy cells approach                                    |
c-----------------------------------------------------------------------
      integer natmx,J,L,NITR,N,K
      double precision GRIDX,GRIDY,GRIDZ,ATWT,AIJ(natmx,natmx),
     & VxIJ(natmx,natmx),VyIJ(natmx,natmx),VzIJ(natmx,natmx),Xp,Yp,Zp,
     & RIJ(natmx,natmx),QQQ,Rp(natmx),WWT(natmx),cutoffIJ,cutoffJI,muij,
     & nuij,Wattot

c---------Step 1 - definition of the position of grid point ------------
c    with rescpect to atomic center J   --------------------------------

          DO J=1,natmx
             Xp = VxIJ(K,J) + GRIDX
             Yp = VyIJ(K,J) + GRIDY
             Zp = VzIJ(K,J) + GRIDZ
             QQQ = Xp**2.d0 + Yp**2.d0 + Zp**2.d0
             Rp(J) = DSQRT(QQQ)
             WWT(J) = 1.d0
          enddo

c----------Step 2 - definition of cut-off function-----------------------
          DO J=1,natmx

c   -------Check if atom is real------- 
          IF(DABS(AIJ(J,1)+1.d0).lt.1d-5) then
          WWT(J) = 0.d0 
          goto 1234
          endif 
c   ------------------------------------
             DO L=1,J-1
              if(ABS(AIJ(K,J)-1.d0).lt.1d-5) then
              cutoffIJ = 1.d0
              goto 1233
              endif
                muij = ( Rp(J) - Rp(L) ) / RIJ(J,L)
                nuij = muij + AIJ(J,L)*(1.d0 - muij**2.d0)

c---------Here nuij - is a generaralized border function for pair ij----------
c----Generation of cutoff function for ij pair with fuzzy polynomial of Becke |
c------------------------------------------------------------------------------

               QQQ = nuij
                DO N=1,NITR
                 QQQ = QQQ*(1.5d0 - 0.5d0*QQQ*QQQ)                 
                ENDDO
                cutoffIJ = 0.5d0*(1.d0-QQQ)
                cutoffJI = 0.5d0*(1.d0+QQQ)
                WWT(J) = cutoffIJ*WWT(J)
                WWT(L) = cutoffJI*WWT(L) 

1233         continue
             ENDDO
1234      continue
          enddo

             Wattot = 0.d0
             DO J=1,natmx
             Wattot = Wattot + WWT(J) 
             enddo

             ATWT = WWT(K)/Wattot
c             write(6,*) 'ATWT =',ATWT,'nnnnnn'
      return
      end 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine make_aij(AIJ,VxIJ,VyIJ,VzIJ,RIJ,atchar,atcoord,bsrd,
     & natmx)
      implicit none
c-----------------------------------------------------------------------
c     This subroutine construct a_ij array, which is necessary         |
c     for construction of Becke fuzzy cells and interatomic            | 
c            pair vectors  V_ij                                        |
c-----------------------------------------------------------------------
      integer natmx,I,J,atchar(natmx),a1,a2
      double precision AIJ(natmx,natmx),bsrd(105),chi,uij,
     & VxIJ(natmx,natmx),VyIJ(natmx,natmx),VzIJ(natmx,natmx),
     & atcoord(3,natmx),RIJ(natmx,natmx),QQQ

        DO I=1,natmx
         a1 = atchar(I)
         AIJ(I,I) = 0.d0
          IF(a1.eq.0) then
          AIJ(I,1) = -1.d0
          goto 4800
          endif 
         
          DO J=1,natmx
          IF(I.EQ.J) goto 4700
           a2 = atchar(J)
           IF(a2.eq.0) THEN
           AIJ(I,J) = 1.d0
           goto 4700 
           ENDIF

           chi = bsrd(a1)/bsrd(a2)
           uij = (chi - 1.d0)/(chi + 1.d0 )
           AIJ(I,J) = uij/(uij*uij-1.d0)

           if(AIJ(I,J).gt.0.5d0) AIJ(I,J)=0.5d0
           if(AIJ(I,J).lt.-0.5d0) AIJ(I,J)=-0.5d0 

4700      continue
          ENDDO
4800      continue
        ENDDO

c-------Interatomic pair vectors----------------------------------------
        DO I=1,natmx
        VxIJ(I,I) = 0.d0
        VyIJ(I,I) = 0.d0
        VzIJ(I,I) = 0.d0        
        RIJ(I,I) = 0.d0
          DO J=1,I-1
          VxIJ(I,J) = atcoord(1,I) - atcoord(1,J) 
          VyIJ(I,J) = atcoord(2,I) - atcoord(2,J)
          VzIJ(I,J) = atcoord(3,I) - atcoord(3,J)
          VxIJ(J,I) = - VxIJ(I,J)
          VyIJ(J,I) = - VyIJ(I,J)
          VzIJ(J,I) = - VzIJ(I,J)
          QQQ = VxIJ(I,J)**2.d0 + VyIJ(I,J)**2.d0 + VzIJ(I,J)**2.d0
          RIJ(I,J) = DSQRT(QQQ)
          RIJ(J,I) = RIJ(I,J)  
          ENDDO 
        ENDDO
      return
      end 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine procgrid(GRIDX,GRIDY,GRIDZ,GRIDW,
     & BLSIZE,ind,ind1,ind2,funit,CUTLIM)

      implicit none

      include 'F12_gr.h'
c-----------------------------------------------------------------------
c       This subroutine copy points and weights to auxiliary array     |
c         and dump it to disk                                          |
c-----------------------------------------------------------------------
      integer BLSIZE,ind,funit,K,ind1,ind2
      double precision GRIDX,GRIDY,GRIDZ,GRIDW,CUTLIM 
     
        ind2 = 0        

        IF(GRIDW.ge.CUTLIM) then
        ind = ind + 1
        ind1 = ind1 + 1 
      PTX(ind1) = GRIDX
      PTY(ind1) = GRIDY
      PTZ(ind1) = GRIDZ
      PTW(ind1) = GRIDW       
c        write(6,*) PTX(ind1), PTY(ind1), PTZ(ind1), '<++==',ind1 
        endif       

c           IF(ind.eq.BLSIZE) then
c           CALL Gwrite(funit,PTX,BLSIZE)
c           CALL Gwrite(funit,PTY,BLSIZE)
c           CALL Gwrite(funit,PTZ,BLSIZE)
c           CALL Gwrite(funit,PTW,BLSIZE) 


c           ind = 0
c           ELSE
c           ind2 = 1
c           ENDIF 

      
      return
      end 

c----------------------------------------------------------------------

c----------------------------------------------------------------------
c-------------------Auxiliary functions--------------------------------

      subroutine Gread(N,ARR,LEN)
      integer N,LEN
      double precision ARR(LEN)
c--------read array ARR of size LEN from the file of unit N ------------
            READ(N) ARR
      return
      end

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine Gwrite(N,ARR,LEN)
      integer N,LEN
      double precision ARR(LEN)
c--------read array ARR of size LEN from the file of unit N ------------
            WRITE(N) ARR
      return
      end


c------------------------------------------------------------------------


 
c--------some auxiliary functions-------------------------------------

      function dref (n, a)
      integer n
      double precision a (*), dref
      dref = a (n)
      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      function iref (n, i)
      integer n, i (*), iref
      iref = i (n)
      return
      end
c---------------------------------------------------------------------

c------------------------------------------------------------------------
       subroutine recharge(natmx,atchar)
       implicit none   
       integer I,natmx,atchar(natmx)

          do I=1,natmx
           IF(atchar(I).eq.110) atchar(I)=1
          enddo
       return
       end 
