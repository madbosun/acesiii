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
      subroutine print_symbolic_constants()
      
      implicit none
      include 'symbolic_constants.h'

      print *, 'norb = ', symbolic_constant_table(1)
      print *, 'nocc = ', symbolic_constant_table(2)
      print *, 'nvirt = ', symbolic_constant_table(3)
      print *, 'bocc = ', symbolic_constant_table(4)
      print *, 'eocc = ', symbolic_constant_table(5)
      print *, 'bvirt = ', symbolic_constant_table(6)
      print *, 'evirt = ', symbolic_constant_table(7)
      print *, 'naocc = ', symbolic_constant_table(8)
      print *, 'nbocc = ', symbolic_constant_table(9)
      print *, 'navirt = ', symbolic_constant_table(10)
      print *, 'nbvirt = ', symbolic_constant_table(11)
      print *, 'baocc = ', symbolic_constant_table(12)
      print *, 'bbocc = ', symbolic_constant_table(13)
      print *, 'eaocc = ', symbolic_constant_table(14)
      print *, 'ebocc = ', symbolic_constant_table(15)
      print *, 'bavirt = ', symbolic_constant_table(16)
      print *, 'bbvirt = ', symbolic_constant_table(17)
      print *, 'eavirt = ', symbolic_constant_table(18)
      print *, 'ebvirt = ', symbolic_constant_table(19)
      print *, 'noccorb = ', symbolic_constant_table(20)
      print *, 'nvirtorb = ', symbolic_constant_table(21)
      print *, 'boccorb = ', symbolic_constant_table(22)
      print *, 'eoccorb = ', symbolic_constant_table(23)
      print *, 'bvirtorb = ', symbolic_constant_table(24)
      print *, 'evirtorb = ', symbolic_constant_table(25)
      print *, 'naoccorb = ', symbolic_constant_table(26)
      print *, 'nboccorb = ', symbolic_constant_table(27)
      print *, 'navirtorb = ', symbolic_constant_table(28)
      print *, 'nbvirtorb = ', symbolic_constant_table(29)
      print *, 'baoccorb = ', symbolic_constant_table(30)
      print *, 'bboccorb = ', symbolic_constant_table(31)
      print *, 'eaoccorb = ', symbolic_constant_table(32)
      print *, 'eboccorb = ', symbolic_constant_table(33)
      print *, 'bavirtorb = ', symbolic_constant_table(34)
      print *, 'bbvirtorb = ', symbolic_constant_table(35)
      print *, 'eavirtorb = ', symbolic_constant_table(36)
      print *, 'ebvirtorb = ', symbolic_constant_table(37)
      print *, 'cc_iter = ', symbolic_constant_table(38)
      print *, 'cc_hist = ', symbolic_constant_table(39)
      print *, 'cc_beg = ',  symbolic_constant_table(40)
      print *, 'scf_iter = ',   symbolic_constant_table(41)
      print *, 'scf_hist = ',  symbolic_constant_table(42)
      print *, 'scf_beg = ',  symbolic_constant_table(43)
      print *, 'ncenters = ', symbolic_constant_table(44)
      print *, 'itrips = ',  symbolic_constant_table(45)
      print *, 'itripe = ',  symbolic_constant_table(46)
      print *, 'ihess1 = ',  symbolic_constant_table(47)
      print *, 'ihess2 = ', symbolic_constant_table(48)
      print *, 'jhess1 = ',  symbolic_constant_table(49)
      print *, 'jhess2 = ',  symbolic_constant_table(50)
      print *, 'subb = ', symbolic_constant_table(51)
      print *, 'sube = ', symbolic_constant_table(52)
      print *, 'sip_sub_segsize = ', symbolic_constant_table(53)
      print *, 'sip_sub_occ_segsize = ',  symbolic_constant_table(54)
      print *, 'sip_sub_virt_segsize ', symbolic_constant_table(55)
      print *, 'sip_sub_ao_segsize = ', symbolic_constant_table(56)

      return 
      end
