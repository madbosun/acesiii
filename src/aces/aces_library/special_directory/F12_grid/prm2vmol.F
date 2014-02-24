      subroutine prm2vmol (
     &     natm, nprm, nzao, nbas, maxshl,
     &     nangatm, nzaoshl, cmpzmt,
     &     prmzao, scrpxz, prm2ao)
c
c This routine transforms primitive coefficients for AO's in ZMAT 
c representation (with a certain order of high-order orbital momentum
c AO's) that are obtained through LIBRA3/A3PRIM to primitive coefficients
c for AO's in "symmetry-adapted" ("VMOL", "computational") representation.
c
c To do this we first sort coefficients for the high-order orbital momentum
c orbitals (A3PRIM gives them in "p_x p_x p_x p_y p_y p_y ...", the SA order
c is c "p_x p_y p_z p_x p_y p_z"), and then transform to the SA order using
c ZMAT2CMP, which is build in V2JA
c
c Igor Schweigert, 2003
c $Id: prm2vmol.F,v 1.2 2004/12/19 20:26:33 igor Exp $
c
      implicit none
c
c The input parameters 
c
      integer
     &     natm, nprm, nzao, nbas, maxshl, nangatm (natm), iintfp,
     &     nzaoshl (maxshl, natm),i,j
c
      double precision
     &     cmpzmt (nzao, nbas), prmzao (nprm, nzao), 
     &     scrpxz (nprm, nzao), prm2ao (nprm, nbas)
c
c Local variables
c
      integer
     &     iatm, iang, nshl, ishl, nza1, izao,
     &     shlang, vaoang, zaooff
c
c Convert the higher angular momentum primitive vectors
c from ZMAT ordering (e.g. ssssxxxyyyzzz...) to VMOL ("computational", 
c "symmetry adapted") ordering (e.g. ssssxyzxyzxyz...)
c
c Variables used in the transformation:
c    ZMAT: s s s s | x x x y y y z z z | ...
c    iang:   1     |         2         | 3
c  shlang:   1     |   2     3     4   | 5
c    nshl:   1     |         3         | 6
c    ishl:   1     |   1     2     3   | ...
c    nza1:   4     |   3     3     3   | ... 
c    izao: 1 2 3 4 | 1 2 3 1 2 3 1 2 3 | ...
c  zaooff: 0       | 4     7     10    | 13 
c    VMOL: s s s s | x y z x y z x y z | ...
c  vaoang: 0       | 4                 | 13
c       write(6,*) 'nprm=',nprm,'nzao=',nzao
c
      vaoang = 0
      zaooff = 0
      do iatm=1, natm
         shlang = 1
         do iang=1, nangatm (iatm)
            nza1 = nzaoshl (shlang, iatm)
            nshl = iang*(iang+1)/2
            do ishl=1, nshl
               do izao=1, nza1
                  call dcopy (
     &                 nprm,
     &                 prmzao (1, zaooff + izao), 1,
     &                 scrpxz (1, vaoang + (izao-1)*nshl + ishl), 1)
               enddo
               zaooff = zaooff + nza1
            enddo
            shlang = shlang + nshl
            vaoang = vaoang + nza1 * nshl 
         enddo
      enddo
c
c Transform reordered primitive coefficients for ZMAT AO's to
c coefficients for VMOL AO's (see VMOL2JA/V2JA.F for details 
c about ZMAT2CMP)
c-------------------------------------------------------------

c      call zero(cmpzmt,nbas*nzao)

c      do i=1,nbas
c      cmpzmt(i,i) = 1.d0
c      enddo


      call xgemm (
     &     'n', 'n', nprm, nbas, nzao,
     &     1.d0, scrpxz, nprm,
     &     cmpzmt, nzao,
     &     0.d0, prm2ao, nprm)
c


      return
      end
