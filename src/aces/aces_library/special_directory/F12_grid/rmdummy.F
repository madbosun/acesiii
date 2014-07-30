      subroutine rmdummy (
     &     natmx, natm, compnorb, comppopv, compmemb, 
     &     atmchr, atmmas, atmcrd)
c
c This routine removes dummy (but not ghost!) atoms from the arrays that
c contain information about the molecule.
c
c Note that all molecular info arrays written into JOBARC by JODA
c include dummy and ghost atoms.
c
c Similar routine are LIBRA3/A3REMOVEDUMMY and HFDFT/REMOVE.
c
c Note that the proper way to go through all the atoms is to loop over
c 1..COMPNORB and 1..COMPPOPV, i.e. through the symmetry ordered atoms,
c since VMOL has already taken care of dummy atoms while building these
c arrays. This routine is needed only if later on somebody decide to
c loop over atoms just by cycling over 1..NATM, thus including the dummy
c atoms.
c
c Igor Schweigert, Feb 2004
c $Id: rmdummy.F,v 1.1 2004/03/01 18:38:30 igor Exp $
c 
      implicit none
c
c     Arguments
c     
      integer
     &     natmx, natm, compnorb, comppopv (compnorb),
     &     compmemb (natmx), atmchr (natmx)
c     
      double precision
     &     atmmas (natmx), atmcrd (3, natmx)
c
c     Local variables
c     
      integer i, j, k, l, n, ndum
c
c     Messages to the user
c     
 5900 format (
     &     t3,'@RMDUMMY-F: The JOBARC record NREALATM does not ',/
     &     t3,'correspond to the number of non-dummy atoms.',/)
c
c     Remove dummy (but not ghost!) atoms from the arrays. Note that
c     COMPPOPV does not count dummy atoms already, but COMPMEMB points
c     to the arrays of atoms including dummy atoms. So COMPMEMB has to
c     be modified as well. Note that we do not modify FULLMEMB.
c     
      j = 0
      ndum = 0
      do i = 1, natmx
         if (atmchr (i) .ne. 0) then
            j = j + 1
            atmmas (j) = atmmas (i)
            atmchr (j) = atmchr (i)
            call dcopy (3, atmcrd (1, i), 1, atmcrd (1, j), 1)
         else
            ndum = ndum + 1
            n = 1
            do k = 1, compnorb
               do l = 1, comppopv (k)
                  if (compmemb (n) + ndum .gt. i)
     &                 compmemb (n) = compmemb (n) - 1
                  n = n + 1
               enddo
            enddo
         endif
      enddo
c
      if (j.ne.natm) then
         write (6,5900)
         call errex
      endif
c    
c          do i=1,natmx
c          write(6,*) atmcrd (1, i), atmcrd (2, i), atmcrd (3, i)
c          enddo
      return
      end
