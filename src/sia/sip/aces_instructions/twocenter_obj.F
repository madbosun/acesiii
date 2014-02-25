      subroutine twocenter_obj(nbas,scfevec,g,a1,a2,MOVAL,obj_type)
      IMPLICIT NONE
c--------------------------------------------------------------------|
c This subroutine calculates values of molecular orbitals and its    |
c derivatives at current grid point g                                |
c--------------------------------------------------------------------|

      include 'mpif.h'
      include 'machine_types.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'dbugcom.h'
      include 'parallel_info.h'
      include 'F12.h'
      include 'F12_gr.h'
      include 'sial_config_params.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

c-----------Variables------------------------------------------------

c--------------Input parameters--------------------------------------
       integer nbas,obj_type,a1,a2,g
       double precision scfevec(nbas,nbas)
c--------------Output array------------------------------------------
       double precision MOVAL(a1:a2)  
c--------------------------------------------------------------------

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff),
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers),
     * zn(max_centers),CONTRVAL(1000),CONTRVALGRA(1000,3),AOVAL(1000),
     * AOVALGRA(1000,3)
      double precision GPRIMVAL(1000),GPRIMVALGRA(1000,3)
      integer angfct(3,100)
c------------------------------------------------------------------------
      integer*8 itop,itop1,itop2,itop3
      integer NEED,NEED1,ierr,iq,m,n,nalpha,ncoeff,angcoef,nprmtot,
     & offset,NXYZ,NRY,NROWMX,NAOShl,i,j,k,l,mm
      double precision gX,gY,gZ,atom_x,atom_y,atom_z
C------------------------------------------------------------------------
      double precision dbuf(1),dref
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


c----------set up angular prefactors--------------------------
           call setangfct(100, angfct)          

c---------Step 1 - get coordinates of current grid point----------------------------------
           gX  = PTX(g)
           gY  = PTY(g)
           gZ  = PTZ(g)
 
C----------Loop over shells --------------------------------------------------------------
           offset = 0
           call zero(MOVAL,a2-a1+1)
                  DO m =1, nshells
C----------------get coordinates of atom of current shell--------------------------------- 
           atom_x = dbuf(icoord+3*(m-1))
           atom_y = dbuf(icoord+3*(m-1)+1)
           atom_z = dbuf(icoord+3*(m-1)+2)

c----------------prepare exponential parameters and contraction coefficients--------------

            call pack_coeffs_2co(dbuf(ialpha), ixalpha, dbuf(ipcoeff),
     *                       ixpcoef, ncfps, npfps, m,
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ibuf(iccbeg), ibuf(iccend),
     *                       indx_cc,
     *                       ccbeg_pack, ccend_pack,
     *                       max_dim_coeff)

c----Total number of primitives (with all angular components) on current shell---------
             angcoef = (ivAngMom(m) + 1)*(ivAngMom(m) + 2)/2
             nprmtot = nPFpS(m)*angcoef
             NAOShl = nCFpS(m)*(2*ivAngMom(m) + 1)
c----Calculate values of primitive gaussians or its gradient at current grid point-----
             call zero(GPRIMVAL,1000)
             call zero(GPRIMVALGRA,3000) 
             call zero(CONTRVAL,1000)
             call zero(CONTRVALGRA,3000)
             call zero(AOVAL,1000)
             call zero(AOVALGRA,3000)

               if(obj_type.eq.1) then
               call valprmgrid (gX,gY,gZ,atom_x,atom_y,atom_z,
     & ivAngMom(m),angfct,alpha_pack,nalpha,nprmtot,GPRIMVAL) 
               else
               call valprmgragrid (gX,gY,gZ,atom_x,atom_y,atom_z,
     & ivAngMom(m),angfct,alpha_pack,nalpha,nprmtot,GPRIMVALGRA)
               endif

c-------------------------------------------------------------------------------------|
c     Prepare contracted functions: multiply primitives values by contraction         |
c     coefficients and corresponding noramlization factors. Then transform contracted |
c     catesion basis functions into the spherical harmonics                           | 
c-------------------------------------------------------------------------------------|
      NXYZ = angcoef
      NRY = 2*ivAngMom(m) + 1
      NROWMX =  ( (ivAngMom(m)/2 + 1) * (ivAngMom(m)/2 + 2) / 2 )
c--------------------------------------------------------------------------------------


               call prepcontrcoeff(nPFpS(m),nCFpS(m),ivAngMom(m),
     & nprmtot,pcoeff_pack,alpha_pack,GPRIMVAL,GPRIMVALGRA,CONTRVAL,
     & CONTRVALGRA,NXYZ,NRY,NROWMX,AOVAL,AOVALGRA)

c-------------------------------------------------------------------------------------|
c    Contract values of AO with the corresponding MO coefficients in order to get     |
c    values of MO at current grid point g                                             | 
c-------------------------------------------------------------------------------------|


c=====================================================================================|
c--------------------Molecular orbitasl values----------------------------------------|
c=====================================================================================|

            if(obj_type.eq.1) then
               do i=a1,a2
                mm = 0 
c-------------------Special MO ordering for p-shell ----------------------------------|
                       IF(ivAngMom(m).eq.1) then
                        do j = 1,nCFpS(m)
                         MOVAL(i) = MOVAL(i) + 
     &                AOVAL(mm+1)*scfevec(mm+1+offset,i) + 
     &                AOVAL(mm+2)*scfevec(mm+3+offset,i) +
     &                AOVAL(mm+3)*scfevec(mm+2+offset,i)
                        mm = mm + 3 
                        enddo
                        goto 222
                       ENDIF

                do j = 1,NAOShl
                mm = mm +1
                MOVAL(i) = MOVAL(i) + AOVAL(mm)*scfevec(mm+offset,i)
                enddo
222             continue
              enddo
            endif
c======================================================================================


c=====================================================================================|
c-------------------x-derivatives of molecular orbitals-------------------------------|
c=====================================================================================|

            if(obj_type.eq.2) then
               do i=a1,a2
                mm = 0
c-------------------Special MO ordering for p-shell-----------------------------------|
                       IF(ivAngMom(m).eq.1) then
                        do j = 1,nCFpS(m)
                         MOVAL(i) = MOVAL(i) +
     &                AOVALGRA(mm+1,1)*scfevec(mm+1+offset,i) +
     &                AOVALGRA(mm+2,1)*scfevec(mm+3+offset,i) +
     &                AOVALGRA(mm+3,1)*scfevec(mm+2+offset,i)
                        mm = mm + 3
                        enddo
                        goto 333
                       ENDIF

                do j = 1,NAOShl
                mm = mm +1
                MOVAL(i) = MOVAL(i)+AOVALGRA(mm,1)*scfevec(mm+offset,i)
                enddo
333             continue
               enddo
            endif
c=====================================================================================|


c=====================================================================================|
c-------------------y-derivatives of molecular orbitals-------------------------------|
c=====================================================================================|
            if(obj_type.eq.3) then
               do i=a1,a2
                mm = 0
c-------------------Special MO ordering for p-shell-----------------------------------|
                       IF(ivAngMom(m).eq.1) then
                        do j = 1,nCFpS(m)
                         MOVAL(i) = MOVAL(i) +
     &                AOVALGRA(mm+1,2)*scfevec(mm+1+offset,i) +
     &                AOVALGRA(mm+2,2)*scfevec(mm+3+offset,i) +
     &                AOVALGRA(mm+3,2)*scfevec(mm+2+offset,i)
                        mm = mm + 3
                        enddo
                        goto 444
                       ENDIF

                do j = 1,NAOShl
                mm = mm +1
                MOVAL(i) = MOVAL(i)+AOVALGRA(mm,2)*scfevec(mm+offset,i)
                enddo
444             continue
               enddo
            endif
c=====================================================================================|

c=====================================================================================|
c-------------------z-derivatives of molecular orbitals-------------------------------|
c=====================================================================================|
            if(obj_type.eq.4) then
               do i=a1,a2
                mm = 0
c-------------------Special MO ordering for p-shell-----------------------------------|
                       IF(ivAngMom(m).eq.1) then
                        do j = 1,nCFpS(m)
                         MOVAL(i) = MOVAL(i) +
     &                AOVALGRA(mm+1,3)*scfevec(mm+1+offset,i) +
     &                AOVALGRA(mm+2,3)*scfevec(mm+3+offset,i) +
     &                AOVALGRA(mm+3,3)*scfevec(mm+2+offset,i)
                        mm = mm + 3
                        enddo
                        goto 555
                       ENDIF

                do j = 1,NAOShl
                mm = mm +1
                MOVAL(i) = MOVAL(i)+AOVALGRA(mm,3)*scfevec(mm+offset,i)
                enddo
555             continue
               enddo
            endif
c=====================================================================================|







c--------------------------------------------------------------------------------------
c--------------------------------------------------------------------------------------


               offset = offset + NAOShl 
c--------------------------------------------------------------------------------------



c-----------------End of loop over shells---------------------------------------------|
                  ENDDO

      return
      end

c---------------------------------------------------------------------------------------
c--------------------------------------------------------------------------------------- 

      subroutine valprmgrid (gX,gY,gZ,atom_x,atom_y,atom_z,angmom,
     & angfct,alpha_pack,nalpha,nprmtot,GPRIMVAL)
      IMPLICIT NONE
c-------------------------------------------------------------------------------------|
c  This subroutine calculates values of gaussian primitives at current grid point     |
c--------------------------------------------------------------------------------------
      integer i,j,k,l,nalpha,nprmtot,angfct (3,100),angmom
      integer max_dim_coeff,angcoef,iprm,ind,shift(0:19)
      parameter (max_dim_coeff = 5000)
      double precision alpha_pack(max_dim_coeff),gX,gY,gZ,atom_x,atom_y,
     & atom_z 
      double precision dist (3), dst, pre,GPRIMVAL(1000)

         angcoef = (angmom + 1)*(angmom + 2)/2

         shift(0)=0
         do i=0,19
         shift(i+1) = shift(i) + (i+1)*(i + 2)/2 
         enddo

         dist (1) = gX - atom_x
         dist (2) = gY - atom_y
         dist (3) = gZ - atom_z
         dst = dist (1)**2.d0 + dist (2)**2.d0 + dist (3)**2.d0
c---------------------------------------------------------------------------------------
         iprm = 0



         do i=1,nalpha
           
         do j=1, angcoef
         iprm = iprm + 1
         
         ind = shift(angmom)+j

            pre =
     &           dist (1)**angfct (1, ind)*
     &           dist (2)**angfct (2, ind)*
     &           dist (3)**angfct (3, ind)



         GPRIMVAL(iprm) = 
     &                    pre * dexp (-1.d0 * alpha_pack (i) * dst)

         enddo
         enddo 

      return
      end

c--------------------------------------------------------------------------------------------
c--------------------------------------------------------------------------------------------
      subroutine valprmgragrid (gX,gY,gZ,atom_x,atom_y,atom_z,angmom,
     & angfct,alpha_pack,nalpha,nprmtot,GPRIMVALGRA)
      IMPLICIT NONE
c-------------------------------------------------------------------------------------|
c  This subroutine calculates values of gaussian primitives derivatives at            |
c                           current grid point                                        |
c-------------------------------------------------------------------------------------|
      integer i,j,k,l,nalpha,nprmtot,angfct (3,100),angmom
      integer max_dim_coeff,angcoef,iprm,ind,shift(0:19),ixyz
      parameter (max_dim_coeff = 5000)
      double precision alpha_pack(max_dim_coeff),gX,gY,gZ,atom_x,atom_y,
     & atom_z,e
      double precision dist (3), dst, pre0(3),pre1(3),
     &  GPRIMVALGRA(1000,3)
      integer pow0 (3,3), pow1 (3,3)

      data
     &     pow0 /-1,  0,  0,  0, -1,  0,  0,  0, -1/,
     &     pow1 / 1,  0,  0,  0,  1,  0,  0,  0,  1/



         angcoef = (angmom + 1)*(angmom + 2)/2

         shift(0)=0
         do i=0,19
         shift(i+1) = shift(i) + (i+1)*(i + 2)/2
         enddo

         dist (1) = gX - atom_x
         dist (2) = gY - atom_y
         dist (3) = gZ - atom_z
         dst = dist (1)**2.d0 + dist (2)**2.d0 + dist (3)**2.d0
c---------------------------------------------------------------------------------------

         iprm = 0

         do i=1,nalpha
         do j=1, angcoef
         iprm = iprm + 1

         ind = shift(angmom)+j

c---------------Angular prefactor derivatives-------------------------------------------
                do ixyz = 1,3
                pre0 (ixyz) = 0d0
                if (angfct (ixyz,ind) .ne. 0 ) pre0 (ixyz) =
     &              dble(angfct (ixyz,ind)) *
     &              dist (1)**(angfct(1,ind)+pow0(1,ixyz))*
     &              dist (2)**(angfct(2,ind)+pow0(2,ixyz))*
     &              dist (3)**(angfct(3,ind)+pow0(3,ixyz))
               pre1 (ixyz) = -2d0 *
     &              dist (1)**(angfct(1,ind)+pow1(1,ixyz))*
     &              dist (2)**(angfct(2,ind)+pow1(2,ixyz))*
     &              dist (3)**(angfct(3,ind)+pow1(3,ixyz))
                enddo
c--------------------------------------------------------------------------------------
                e = dexp ( - alpha_pack (i) * dst) 
               
                do ixyz=1,3
                 
         GPRIMVALGRA(iprm,ixyz) =
     &   e * ( pre0 (ixyz) + pre1 (ixyz) * alpha_pack (i) )
              
                enddo

         enddo
         enddo

      return
      end

c---------------------------------------------------------------------------------------
c---------------------------------------------------------------------------------------
      subroutine prepcontrcoeff(nprm,ncontr,angmom,nprmtot,pcoeff_pack,
     &alpha_pack,GPRIMVAL,GPRIMVALGRA,CONTRVAL,CONTRVALGRA,NXYZ,NRY,
     & NROWMX,AOVAL,AOVALGRA)
      IMPLICIT NONE
c-------------------------------------------------------------------------------------|
c This subroutine prepares contracted basis functions: multiply them by noramlization |
c          factors  and transforms to basis of spherical harmonics                    |
c-------------------------------------------------------------------------------------|
      integer nprm,ncontr,i,j,k,l,dft,angcoef,mm,ind,NXYZ,NRY,NROWMX 
      integer max_dim_coeff,angmom,nprmtot,NROW(NRY)
      parameter (max_dim_coeff = 5000)
      double precision pcoeff_pack(max_dim_coeff),NormFactor(nprm),
     & alpha_pack(max_dim_coeff),contrcoeff(ncontr,nprmtot),
     & TEMP(NXYZ)
      integer  ROW(NROWMX,NRY)
      double precision TMAT( NROWMX,NRY)
      double precision pi,NOM,DEN,qq,qq1,tt,ppx,ppy,ppz,GPRIMVAL(1000),
     & GPRIMVALGRA(1000,3),CONTRVAL(1000),CONTRVALGRA(1000,3),
     & AOVAL(1000),AOVALGRA(1000,3)
      parameter (pi = 3.141592653589793D0)         
c--------------Generation of normalization factors-------------------------------------
            DEN = DBLE(dft(2*angmom+1)) * DSQRT(pi)
            qq = 2.d0**( 2.d0 * DBLE(angmom) + 3.d0 + 0.5d0)
            qq1 = (2.d0 * DBLE(angmom) + 3.d0)/2.d0
              do i=1,nprm     
              NOM = qq * alpha_pack(i)**qq1
              NormFactor(i) = DSQRT(NOM/DEN)
              enddo
c--------------------------------------------------------------------------------------
            angcoef = (angmom + 1)*(angmom +2) / 2
            k=0
            do i=1,ncontr
            mm = 0  
            do j=1,nprm
            k=k+1
            pcoeff_pack(k) = pcoeff_pack(k)*NormFactor(j)/dsqrt(pi)
                do l=1,angcoef
                contrcoeff(i,mm+l) = pcoeff_pack(k) 
                enddo
            mm = mm + angcoef
            enddo
            enddo

c----Calculation of the values of contracted functions and its derivatives-------------

            k=0
            do i=1,ncontr

             do l=1,angcoef
              k=k+1
              tt = 0.d0
              ppx = 0.d0
              ppy = 0.d0
              ppz = 0.d0   
              ind = 0
c--------------------------------                         
              do j=1,nprm
              tt = tt + GPRIMVAL(ind+l)*contrcoeff(i,ind+l)
              ppx = ppx + GPRIMVALGRA(ind+l,1)*contrcoeff(i,ind+l)
              ppy = ppy + GPRIMVALGRA(ind+l,2)*contrcoeff(i,ind+l)
              ppz = ppz + GPRIMVALGRA(ind+l,3)*contrcoeff(i,ind+l)
              ind = ind + angcoef
              enddo
c----------------------------------
              CONTRVAL(k) = tt
              CONTRVALGRA(k,1) = ppx
              CONTRVALGRA(k,2) = ppy
              CONTRVALGRA(k,3) = ppz
             enddo

            enddo

c#######################################################################################
c-----------------Transformation of cartesian AO to sherical AO-------------------------
c#######################################################################################


c--------Step 1 - generation of the transformation matrix-------------------------------

       CALL zero(ROW,NRY*NROWMX)
       CALL zero(TMAT,NRY*NROWMX)


       call  cart2sphtran( NXYZ,NRY,NROWMX,angmom,TEMP, 
     +
     +                            NROW,
     +                            ROW,
     +                            TMAT )

c--------Step 2 - transformation of cartesian AO to spherical  AO----------------------

            k=0
            mm = 0 
            do i=1,ncontr
             do l=1,NRY
              k=k+1
              tt = 0.d0
              ppx = 0.d0
              ppy = 0.d0
              ppz = 0.d0


                do j=1,NROWMX
                     ind = ROW(j,l)
                     IF(ind.ne.0) then 
                tt = tt + TMAT(j,l)*CONTRVAL(mm+ind)
                ppx = ppx + TMAT(j,l)*CONTRVALGRA(mm+ind,1)     
                ppy = ppy + TMAT(j,l)*CONTRVALGRA(mm+ind,2)
                ppz = ppz + TMAT(j,l)*CONTRVALGRA(mm+ind,3)
                     ENDIF
                enddo
c----------------------------------
              AOVAL(k) = tt
              AOVALGRA(k,1) = ppx
              AOVALGRA(k,2) = ppy
              AOVALGRA(k,3) = ppz
             
             enddo
             mm = mm + angcoef
           enddo

c--------------------------------------------------------------------------------------

      return
      end

c---------------------------------------------------------------------------------------
c---------------------------------------------------------------------------------------
      function dft(n)
      integer dft,n,k 
c-------------------------------------------------------------|
c This function generates double factorial of integr number n |
c-------------------------------------------------------------|
      dft =1       
      do k=n,1,-2
      dft = dft*k
      enddo
      return
      end

