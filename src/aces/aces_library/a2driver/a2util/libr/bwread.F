
c This routine reads the bwdata file. If this file does not exist, then the
c usual CC is performed.

      subroutine bwread(invcc,iuhf)
      logical invcc
      logical bIAmOne
      character*512 line
      character*80 fname
      real*8 s
      integer ilength,nbas
      character*1 z
#include "bwcc.com"

      bIAmOne = .true.

c initialize warnings
      do i = 1, maxbwwarnings
         bwwarning(i) = .false.
      end do

      intruder = 1.d300
      open(unit=1,file='intruder',
     &     form='formatted',status='old',err=991)
      read(1,*) intruder
      close(1)
      intruder = abs(intruder)
991   continue
      if (intruder.lt.1.d300.and.bIAmOne) then
         write(6,*)'Intruder limit set to',intruder
      end if

      correctiontype = 0
      open(unit=1,file='correctiontype',
     &     form='formatted',status='old',err=992)
      read(1,*) correctiontype
      close(1)
992   continue

      ihubaccorr = 1
      open(unit=1,file='nohubaccorr',status='old',err=994)
      ihubaccorr = 0
      close(1)
994   continue

      ihomotop = 0
      open(unit=1,file='homotop',
     &     form='formatted',access='sequential',status='old',err=993)
      read(1,*) ihomotop
      read(1,*) lambdahomotop
      read(1,*) diishonset
      hfaktmax = 1.d-6
      read(1,*,end=9931,err=9931) hfaktmax
9931  continue
      close(1)
993   continue
      if (ihomotop.gt.0) then
         ihubaccorr = 0
         if (invcc.and.bIAmOne) then
            write(6,*)
            write(6,*) 'homotopic solution attempted, parameters ',
     &                 ihomotop,lambdahomotop,diishonset,hfaktmax
            write(6,*) 'correction type: ',correctiontype
            write(6,*)
         end if
      end if

      bwgossip = .false.
      open(unit=1,file='bwgossip',status='old',err=995)
      bwgossip = .true.
      close(1)
995   continue
      if (bIAmOne) write(6,*)
      if (bIAmOne) write(6,*)

      isbwcc = .true.
      open(unit=1,file='masik',
     &     form='formatted',access='sequential',status='old',err=998)
      if (bIAmOne) then
         write(6,*) 'performing integral processing for masik program'
      end if
      masik = .true.
      close(1)
      goto 111
998   masik = .false.
111   continue

      open(unit=1,file='bwdata',
     &     form='formatted',access='sequential',status='old',err=999)
      if (invcc) then
         if (bIAmOne) write(6,*) 'running mr-bwcc program'
      else
         if (bIAmOne) then
            write(6,*) 'performing integral processing for mr-bwcc pgm'
         end if
      end if

c we need IUHF before crapsi is called
c vscf will create this file
      open(unit=2,file='IUHF',
     &     form='unformatted',access='sequential',status='old')
      read(2) iuhf
      close(2)

      do j = 1, maxref
      do i = 1, maxorb
         iocc(i,j,1) = 0
         iocc(i,j,2) = 0
      end do
      end do
      do i = 1, maxorb
         iocc0(i,1) = 0
         iocc0(i,2) = 0
      end do
      scfrefread = .false.

c here we read the data
      nref = 0
1     read(1,10,end=2) line
10    format(a512)
      if (bIAmOne) write(6,*) 'input reference: ',line
      if (line(1:3).eq.'SCF' .or. line(1:3).eq.'scf') then
c what we have read in previous was just the SCF configuration
c not a reference configuration
         scfrefread = .true.
         do i = 1, maxorb
            iocc0(i,1)     = iocc(i,nref,1)
            iocc0(i,2)     = iocc(i,nref,2)
            iocc(i,nref,1) = 0
            iocc(i,nref,2) = 0
         end do
      else
         nref = nref+1
      end if
      if (nref.gt.maxref) then
         write(*,*) 'too many references, increase MAXREF'
         call aces_exit(1)
      end if

c process the string into array
      do i = 1, maxorb
         z = line(i:i)
         if (z.eq.'2') then
            iocc(i,nref,1) = 1
            iocc(i,nref,2) = 1
         else
            if (z.eq.'a' .or. z.eq.'b') then
               if (iuhf.eq.0) then
                  write(*,*) 'open shell references not allowed for ',
     &                       'RHF based calculation'
                  call aces_exit(1)
               else
                  if (z.eq.'a') then
                     iocc(i,nref,1) = 1
                  else
                     iocc(i,nref,2) = 1
                  end if
               end if
            else
               if (z.ne.'0'.and.z.ne.' ') then
                  if (bIAmOne) then
                     write(6,*) 'unexpected character "',z,'" in ref.'
                     write(6,*) line
                  end if
                  write(*,*) 'illegal character in configuration ',
     &                       'specification'
                  call aces_exit(1)
               end if
            end if
         end if
      end do
      goto 1
2     close(1)
      if (bIAmOne) then
         write(6,*)
         write(6,*) 'number of bw references: ',nref
         write(6,*)
      end if
      goto 100
999   isbwcc = .false.
      if (invcc) then
         if (bIAmOne) write(6,*) ' running standard cc program'
      else
         if (bIAmOne) then
            write(6,*) 'performing integral processing for standard CC'
         end if
      end if
      close(1)
      nref = 1
c in cc program read and store vacuum energies and diagonal fock corrections
c in intprc this does not exist yet, of course
100   continue

      if (.not.isbwcc) goto 200
      if (invcc) then
c this loop has to be parallelized since FOCKCD was generated only for
c appropriate irefs by xintprc
         do iref = 1, nref
            if (mod(iref-1,1).ne.1-1) goto 1234
            call gfname('FOCKCD  ',fname,ilength)
            open(unit=99,file=fname(1:ilength),
     &           form='unformatted',status='old')
            read(99) nbas
            do is = 1, 1+iuhf
               do i = 1, nbas
                  read(99) fockcd(i,iref,is)
               end do
            end do
            read(99) hdiagcontr(iref)
c read also indices for assignment of offdiagonal heff elements
            if (iuhf.eq.0) then
               ibot = 3
            else
               ibot = 1
            end if
            do is = ibot, 3
               read(99) internnum(iref,is)
               read(99) (internfrom(j,iref,is),
     &                   j=1,(maxref*(maxref-1)/2))
               read(99) (internto(j,iref,is),
     &                   j=1,(maxref*(maxref-1)/2))
               read(99) (internindex(j,iref,is),
     &                   j=1,(maxref*(maxref-1)/2))
            end do
            if (iuhf.ne.0) then
               do is = 1, 2
                  read(99) internnum1(iref,is)
                  read(99) (internfrom1(j,iref,is),
     &                      j=1,(maxref*(maxref-1)/2))
                  read(99) (internto1(j,iref,is),
     &                      j=1,(maxref*(maxref-1)/2))
                  read(99) (internindex1(j,iref,is),
     &                      j=1,(maxref*(maxref-1)/2))
               end do
            end if
            close(99)
1234     end do
      else
c initialize indices for assignment of offdiagonal heff elements
         if (iuhf.eq.0) then
            ibot = 3
         else
            ibot = 1
         end if
         do is = ibot, 3
            do iref = 1, nref
               internnum(iref,is) = 0
               do j = 1, maxref*(maxref-1)/2
                  internindex(j,iref,is) = 0
                  internfrom(j,iref,is)  = 0
                  internto(j,iref,is)    = 0
               end do
            end do
         end do
         if (iuhf.ne.0) then
            do is = 1, 2
               do iref = 1, nref
                  internnum1(iref,is) = 0
                  do j = 1, maxref*(maxref-1)/2
                     internindex1(j,iref,is) = 0
                     internfrom1(j,iref,is)  = 0
                     internto1(j,iref,is)    = 0
                  end do
               end do
            end do
         end if
      end if
c the following in fact not necessary
200   continue

c read info which roots we want
      nbwstates   = 1
      ibwstate(1) = 1
      cbwstate(1) = 1.0
      if (isbwcc) then
         open(unit=1,file='bwroot',
     &        form='formatted',access='sequential',status='old',err=997)
333      read(1,'(i3,f10.7)',end=3) ibwstate(nbwstates),
     &                              cbwstate(nbwstates)
         if (ibwstate(nbwstates).le.0 .or.
     &       ibwstate(nbwstates).gt.nref) then
            write(6,*) 'illegal state number requested'
            call aces_exit(1)
         end if
         nbwstates = nbwstates+1
         goto 333
3        close(1)
         nbwstates = nbwstates-1
         if (nbwstates.gt.nref) then
            write(6,*) 'too many states to average'
            call aces_exit(1)
         end if
         s = 0.0
         do i = 1, nbwstates
            s = s + cbwstate(i)*cbwstate(i)
         end do
         if (s.eq.0.0) then
            do i = 1, nbwstates
               cbwstate(i) = 1d0/nbwstates
            end do
         end if
997      continue

         if (invcc.and.bIAmOne) then
            write(6,*) 'BW roots requested: ',nbwstates
            do i = 1, nbwstates
               write(6,*) 'BW root ',ibwstate(i), cbwstate(i)
            end do
            write(6,*)
            write(6,*)
         end if

      end if

      iref = 1

      if (.not.isbwcc) ihubaccorr=0
      if (ihubaccorr.gt.0.and.bIAmOne) then
         write(6,*) 'A posteriori correction will be calculated'
         write(6,*) 'Correction type: ',correctiontype
         write(6,*)
      end if

c check if parallelization makes sense
      if (.not.isbwcc.and.1.gt.1) then
         write(6,*)' PARALLELIZATION IMPLEMENTED ONLY FOR MR-BWCC '
         call aces_exit(1)
      end if
      if (mod(nref,1).ne.0.and.bIAmOne) then
         write(6,*) ' PARALLELIZATION EFFICIENCY WILL NOT BE OPTIMAL'
         write(6,*) ' NUMBER OF PROCESSORS SHOULD BE DIVISOR OF NREF'
      end if
      if (1.gt.nref) then
c it has no sense and broadcast would not work either:
         if (bIAmOne) then
            write(6,*)' IT HAS NO SENSE TO RUN ON MORE PROCESSORS THEN'
            write(6,*)' NUMBER OF REFERENCES. AND IT DOES NOT WORK'
         end if
         call aces_exit(1)
      end if

      return
      end

