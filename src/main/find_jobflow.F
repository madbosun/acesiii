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
      subroutine find_jobflow(calc,dropmo,ref,geom_opt,
     *                                       vib, excite, instab, props,
     *                                       grad_calc, density, 
     *                                       dip_calc, dea_calc, 
     *                                       jobflow)
c-------------------------------------------------------------------------
c   Writes the SIAL_PROGRAM parameters for a default set of SIAL programs
c   determined by the parameters calc, dropmo, ref, geom_opt, and vib.
c   The SIAL_PROGRAM parameters are written to ZMAT.AUTO.
c-------------------------------------------------------------------------
      implicit none
      integer calc,dropmo,ref,geom_opt, vib, excite, instab, props,
     &        density, dip_calc, dea_calc, dip_singlet_root,
     &        dip_triplet_root, dea_singlet_root, dea_triplet_root

      integer grad_calc 
      character*80 jobflow
      integer ierr
      integer n, str_trimlen

c---------------------------------------------------------------------------
c   Determine the jobflow required for the combination of parameters.
c---------------------------------------------------------------------------

c      print *,'ref, calc, dropmo, geom_opt, vib ',
c     *   ref, calc, dropmo, geom_opt, vib, excite
      jobflow = 'UNDEFINED'

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'SCF_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'SCF_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'SCF_RHF_GRADIENT'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'SCF_UHF_GRADIENT'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'SCF_RHF_HESSIAN'     ! analytical hessian calc
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 3) then
         if (grad_calc .eq. 2) then
            jobflow = 'SCF_RHF_ENERGY'    ! numerical gradient calc
         else
            jobflow = 'SCF_RHF_GRADIENT'    ! numerical hessian calc
         endif 
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'SCF_UHF_HESSIAN'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'SCF_UHF_ENERGY'
         else
            jobflow = 'SCF_UHF_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'MP2_RHF_HESSIAN'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'MP2_RHF_ENERGY'
         else
            jobflow = 'MP2_RHF_GRADIENT'
         endif
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'MP2_UHF_HESSIAN'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'MP2_UHF_ENERGY'
         else
            jobflow = 'MP2_UHF_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'MP2_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'MP2_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'MP2_RHF_GRADIENT'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'MP2_UHF_GRADIENT'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_RHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_UHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'LCCSD_RHF_ENERGY'
         else
            jobflow = 'LCCSD_RHF_GRADIENT'
         endif
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'LCCSD_UHF_ENERGY'
         else 
            jobflow = 'LCCSD_UHF_GRADIENT'
         endif  
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         jobflow = 'LCCSD_RHF_DROPMO_GRADIENT'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         jobflow = 'LCCSD_UHF_DROPMO_GRADIENT'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_RHF_ENERGY'
         else
            jobflow = 'CCSD_RHF_GRADIENT'
         endif
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_UHF_ENERGY'
         else
            jobflow = 'CCSD_UHF_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_RHF_DROPMO_ENERGY'
         else
            jobflow = 'CCSD_RHF_DROPMO_GRADIENT'
         endif 
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_UHF_DROPMO_GRADIENT'
         else
            jobflow = 'CCSD_UHF_DROPMO_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 22 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_TRIPLES_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 22 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_TRIPLES_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 22 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_DROPMO_TRIPLES_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 22 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_DROPMO_TRIPLES_ENERGY'
      endif

      n = str_trimlen(jobflow)
      if (instab .gt. 0) jobflow((n+1):(n+7)) = '_INSTAB'

      if (excite .eq. 3) then
         if (ref .eq. 0 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0) then
            jobflow = 'EOM_RHF_CCSD_ENERGY'
         endif

         if (ref .eq. 0 .and.
     *       calc .eq. 10 .and.
     *       dropmo .ne. 0) then
            jobflow = 'EOM_RHF_CCSD_DROPMO_ENERGY'
         endif

         if (ref .eq. 1 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0) then
            jobflow = 'EOM_UHF_CCSD_ENERGY'
         endif

         if (ref .eq. 1 .and.
     *       calc .eq. 10 .and.
     *       dropmo .ne. 0) then
            jobflow = 'EOM_UHF_CCSD_DROPMO_ENERGY'
         endif
C
C   Watson
C
         if (ref .eq. 0 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0  .and.
     *       props  .eq. 1)   then
            jobflow = 'EOM_RHF_CCSD_DENS_ENERGY'
         endif

         if (ref .eq. 1 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0  .and.
     *       props   .eq. 1)  then
            jobflow = 'EOM_UHF_CCSD_DENS_ENERGY'
         endif
C
C   Watson
C
      endif
C
C Ajith Perera, 07/2013. The isotropic Hyperfine coupling 
c Constants.
C  
   
      if (ref .eq. 1 .and. Calc .eq. 0 .and. props .eq. 1) 
     &   jobflow = 'HF_AISO'

      if (ref .eq. 1 .and. Calc .eq. 1 .and. props .eq. 1) 
     &   jobflow = 'MBPT2_AISO'

      if (ref .eq. 1 .and. Calc .eq. 10 .and. props .eq. 1 .and.
     &   density .eq. 0) jobflow = 'CCSD_RELAX_DENS_AISO'

      if (ref .eq. 1 .and. Calc .eq. 22 .and. props .eq. 1 .and.
     &   density .eq. 0) jobflow = 'CCSDT_RELAX_DENS_AISO'

      if (ref .eq. 1 .and. Calc .eq. 10 .and. props .eq. 1 .and.
     &   density .eq. 1) jobflow = 'CCSD_RESPONSE_DENS_AISO'

      if (ref .eq. 1 .and. Calc .eq. 10 .and. props .eq. 1 .and.
     &   density .eq. 0 .and. dropmo .ne. 0) 
     &   jobflow = 'CCSD_RELAX_DENS_DROPMO_AISO'
  
C
C EOM-DIP and DEA (singlets and triplets)
C
C Note that currently singlet states only DIP/DEA calculations
C require UHF refrence and it is not necessary and need to be
C changed (require new SIAL codes).


      call Igetrec(-1,'JOBARC','DIPSYMA ',1, dip_singlet_root)
      call Igetrec(-1,'JOBARC','DIPSYMB ',1, dip_triplet_root)

      call Igetrec(-1,'JOBARC','DEASYMA ',1, dea_singlet_root)
      call Igetrec(-1,'JOBARC','DEASYMB ',1, dea_triplet_root)

      if (ref .eq. 1 .and. Calc .eq. 10 .and.dip_calc .eq. 2 .and.
     &    dip_singlet_root .gt. 0) jobflow = 'EOMDIP_CCSD_SINGLET'

      if (ref .eq. 1 .and. Calc .eq. 10 .and. dip_calc .eq. 2 .and.
     &    dip_triplet_root .gt. 0) jobflow = 'EOMDIP_CCSD_TRIPLET'

      if (ref .eq. 1 .and. Calc .eq. 10 .and. dea_calc .eq. 2 .and.
     &    dea_singlet_root .gt. 0) jobflow = 'EOMDEA_CCSD_SINGLET'

      if (ref .eq. 1 .and. Calc .eq. 10 .and. dea_calc .eq. 2 .and.
     &    dea_triplet_root) jobflow = 'EOMDEA_CCSD_TRIPLET'

      return
      end

