      subroutine F12_grid()
      IMPLICIT NONE
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
c------------Definition of the variables----------------------------------

      integer BLSIZE,RADNUM,ANGNUM,ibsrd,ir_weight,ir_point,ierr,
     & lebed_weight,x_lebed,y_lebed,z_lebed,coretop,NEED

      integer*8 itop
C------------------------------------------------------------------------
      double precision dbuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
#else
      common dbuf
#endif
      double precision vnn

#ifdef ALTIX
      dptr = dshptr   ! set dbuf ptr to shared memory.
#endif

c----------------------------------------------------------------------------


      parameter (BLSIZE = 500)
c---------Number of points in radial and angular grid-------------------------------
         RADNUM = 50
         ANGNUM = 194
c----------------------------------------------------------------------------------
         ibsrd = 1
         ir_weight = ibsrd + 105
         ir_point = ir_weight + RADNUM
         lebed_weight = ir_point + RADNUM
         x_lebed = lebed_weight + ANGNUM
         y_lebed = x_lebed + ANGNUM
         z_lebed = y_lebed + ANGNUM
         NEED = z_lebed + ANGNUM
c-----------------Memory allocation for temporal arrays-----------------------------

      call mem_alloc(dbuf, NEED, bytes_per_double, itop, .true., ierr)
      if(ierr.ne.0) then
      write(6,*) 'Not enough memory!!!!'
      return 
      endif
c--------Generation of the grid using Becke fuzzy cells method---------------------

         call beckegrid( RADNUM,ANGNUM,dbuf(itop+ir_point-1),
     &   dbuf(itop+ir_weight-1),dbuf(itop+ibsrd-1),dbuf(itop+x_lebed-1),
     &   dbuf(itop+y_lebed-1),dbuf(itop+z_lebed-1),
     &   dbuf(itop+lebed_weight-1))

c----------------------------------------------------------------------------------



c----------------------------------------------------------------------------------      

      return
      end 
       
