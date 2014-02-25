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
      subroutine fragment_orbitals(fend_nfps, fnshells, atoms) 
      implicit none 
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'frag_def.h'
      include 'dropmo.h'

      integer i, j, k, natoms, watom, nao, otemp    
      integer fnshells 
      integer fend_nfps(fnshells)
      integer atoms(fnshells) 

c     logical frag_calc 
c     integer nfrags, natoms_frag, watom_frag 
c     common /FRAGDEF/frag_calc, nfrags, natoms_frag(100), 
c    *                watom_frag(100,100) 

c     integer fragAO, n_occ, n_aocc, n_bocc, n_virt, n_avirt, n_bvirt
c     common /FORBS/fragAO(100), n_occ(100), n_aocc(100), n_bocc(100), 
c    *              n_virt(100), n_avirt(100), n_bvirt(100) 
      
c----------------------------------------------------------------------- 
c  Determine the number of AO's in each fragment 
c----------------------------------------------------------------------- 

      do j = 1, nfrags 
         fragAO(j) = 0 
      do k = 1, natoms_frag(j) 
      do i = 1, fnshells
         watom = atoms(i) 
         if (i .eq. 1) nao = fend_nfps(i) 
         if (i .gt. 1) nao = fend_nfps(i) - fend_nfps(i-1) 
         if (watom .eq. watom_frag(j,k)) fragAO(j) = fragAO(j) + nao  
      enddo 
      enddo 
      enddo 
      
c----------------------------------------------------------------------- 
c  Determine the number of occupied orbitals in each fragment.  
c----------------------------------------------------------------------- 

      do i = 1, nfrags 
         n_occ(i)   = 0 
         n_aocc(i)  = 0 
         n_bocc(i)  = 0 
         n_virt(i)  = 0 
         n_avirt(i) = 0 
         n_bvirt(i) = 0 
         otemp = 0 
      do j = 1, natoms_frag(i) 
         watom = watom_frag(i,j) 
      do k = 1, max_centers ! should ne natoms 
         if (k .eq. watom) then 
            otemp      = otemp + charge(k) ! charge(k)/2 
c           n_bocc(i)  = n_bocc(i) + otemp  
c           n_aocc(i)  = n_aocc(i) + charge(k) - otemp  
         endif 
      enddo 
      enddo 
         n_bocc(i)  = otemp /2  
         n_aocc(i)  = otemp - n_bocc(i)   
         n_avirt(i) = fragAO(i) - n_aocc(i) 
         n_bvirt(i) = fragAO(i) - n_bocc(i) 
      enddo 

c Check for dropmo calculation 
      write(6,*) ' FDROP :', fdrop_flag 
      if (fdrop_flag) then 
c         ndropmo_a = 5 
c         ndropmo_b = 5 
         write(6,*) ' NDROP :', ndropmo_a, ndropmo_b 
         n_aocc(1) = ndropmo_a  
         n_aocc(2) = nalpha_occupied - ndropmo_a  
         n_bocc(1) = ndropmo_b  
         n_bocc(2) = nbeta_occupied - ndropmo_b  
      endif 

c     go to 11 
      if (me .eq. 0) then 

      write(6,10) 
      write(6,20) nfrags  
      write(6,30)   
      write(6,40)  
      do i = 1, nfrags 
         write(6,50) i, natoms_frag(i)   
      do j = 1, natoms_frag(i) 
         write(6,60) watom_frag(i,j), watom_frag(i,j) 
      enddo 
      enddo 

      do i = 1, nfrags
         write(6,*) ' FRAGMENT = ', i 
         write(6,*) '    AO   = ', i, fragAO(i)  
         write(6,*) '    OCCA = ', i, n_aocc(i) 
         write(6,*) '    OCCB = ', i, n_bocc(i) 
         write(6,*) '   VIRTA = ', i, n_avirt(i) 
         write(6,*) '   VIRTB = ', i, n_bvirt(i) 
      enddo 

10    format('Parameters of Fragment calculation') 
20    format('Number of Fragments',I8) 
30    format('Atoms in Fragments') 
40    format(6x,'Frag',12x,'Natoms') 
50    format(6x,I8,10x,I8) 
60    format(4x,'ATOM LABEL',I8,10x,'ATOM',I8) 

      endif 
11    continue 

      end 



