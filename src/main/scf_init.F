C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine scf_init(vnn)
c-------------------------------------------------------------------------
c   Handles most set-up chores for a run which begins with no JOBARC
c   in existence.  It is assumed that we must start with an SCF to begin
c   such a job.
c-------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'mxatms.par' 

      integer ierr
      common /flags/ iflags
      integer iflags(100)

      integer nproton, icharg, imult, iuhfx, nrhs, ialpex, numel
      integer nallproton, noff  
      integer dummy, i, nocc_present
      logical*8 l8false
      integer ECP, NATOMS, ICRCOR(MxAtms) 
      integer nallalpha_occupied, nallbeta_occupied 
      double precision vnn 
      double precision nuclear_nuclear_repulsion_energy

      double precision dbuf(1)
      common dbuf

      logical is_open
 
      inquire (file='JOBARC', exist=jobarc_exists)
      if (.not. jobarc_exists) then
         print *,'Running initial xjoda...'
         call c_flush_stdout()

         call joda_main()
         print *,'Initial xjoda is complete...'
         call c_flush_stdout()
      endif

c---------------------------------------------------------------------------
c   Initialize ACES routines.
c---------------------------------------------------------------------------
                                                                                
      call aces_init_rte()
      call aces_ja_init()

c----------------------------------------------------------------------------
c   xjoda sometimes destroys the OCCUPYA0, OCCUPYB0 data on JOBARC for 
c   some obscure reason.  To guard against that possibility, we must 
c   check to see if it is available now.  If not, we calculate it from
c   ourselves and write it to JOBARC.
c----------------------------------------------------------------------------

      call igetrec(0, 'JOBARC', 'OCCUPYA0', nocc_present, dummy)
      if (nocc_present .gt. 0) then

c---------------------------------------------------------------------
c   Go ahead and read the occupied values off JOBARC.
c---------------------------------------------------------------------

         call igetrec(1, 'JOBARC', 'OCCUPYA0', 1, nalpha_occupied)
         call igetrec(1, 'JOBARC', 'OCCUPYB0', 1, nbeta_occupied)

         return
      endif

c---------------------------------------------------------------------------
c   OCCUPYA0/OCCUPYB0 will not exist until after the SCF.  So we have to 
c   calculate the number of occupied alpha and beta spin electrons.
c---------------------------------------------------------------------------

      CALL IGETREC(20,'JOBARC','IFLAGS  ', 100, iflags)
      CALL IGETREC(20,'JOBARC','NMPROTON',1,NPROTON)

C--------------------------------------------------------------------------
c   If this is an ECP calculation then we must modify NPROTON to obtain
c   correct number of occupied orbitals. Why is this not done earlier??   
C--------------------------------------------------------------------------

      ECP = iflags(71) 
      IF (ECP .EQ. 1) THEN 
         NallPROTON = NPROTON 

         CALL SET_ECP_CHARGES(MAX_CENTERS,NATOMS,ICRCOR) 
         CALL IGETREC(20,'JOBARC','ATOMCHRG',NATOMS,ICRCOR)
c 
c        do i = 1, natoms 
c           write(6,*) 'ATOMCHRG  ', i, icrcor(i) 
c        enddo 

c        natoms = 2
c        write(6,*) ' Natoms = ', natoms 

         nproton = 0 
         do i = 1, natoms 
#ifdef _DEBUG_LVL0
            write(6,*) '  ', i, icrcor(i) 
#endif
            nproton = nproton + icrcor(i) 
         enddo 

c        icrcor(1) = 7 
c        icrcor(2) = 7 

         do i = 1, natoms 
            charge(i) = icrcor(i) 
         enddo 

c        nproton = 14 

         noff = NallPROTON - NPROTON 
#ifdef _DEBUG_LVL0
         write(6,*) ' Nproton = ', nallproton, nproton, noff  
#endif 
      ELSE 
      ENDIF 

      ICHARG=IFLAGS(28)
      IMULT=IFLAGS(29)

C--------------------------------------------------------------------------
c   Determine alpha and beta occupied from charge and multiplicity.
C--------------------------------------------------------------------------

      NUMEL=NPROTON-ICHARG
      IALPEX=IMULT-1
      NRHS=NUMEL-IALPEX
      IF(MOD(NRHS,2).NE.0)THEN
         print *,'Error: Specified charge and multiplicity are ',
     &           'impossible.   Try again.'
         call abort_job()
      ENDIF
      nbeta_occupied=NRHS/2
      nalpha_occupied = nbeta_occupied + ialpex
      if (nalpha_occupied .eq. 0) call abort_job()

c     IF (ECP .EQ. 1) THEN 

c        NUMEL=NallPROTON-ICHARG
c        IALPEX=IMULT-1
c        NRHS=NUMEL-IALPEX
c        IF(MOD(NRHS,2).NE.0)THEN
c           print *,'Error: Specified charge and multiplicity are ',
c    &              'impossible.   Try again.'
c           call abort_job()
c        ENDIF
c        nallbeta_occupied=NRHS/2
c        nallalpha_occupied = nallbeta_occupied + ialpex
c        if (nallalpha_occupied .eq. 0) call abort_job()

c     ENDIF 
      
c-------------------------------------------------------------------------
c   Write the values to JOBARC.
c-------------------------------------------------------------------------

      call iputrec(1, 'JOBARC', 'OCCUPYA0', 1, nalpha_occupied)
      call iputrec(1, 'JOBARC', 'OCCUPYB0', 1, nbeta_occupied)

      return
      end
