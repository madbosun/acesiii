
c This routine finalizes the job archive environment (via the JAINDX file)
c in order for getrec and putrec to work properly for other member executables.







      subroutine aces_ja_fin
      implicit none

c INTERNAL VARIABLES
      character*80 szJAINDX
      integer       iJAINDX
      integer i, iStat
      logical bOpened

c COMMON BLOCKS
c jobarc.com : begin

c This data tracks the contents of the JOBARC file. 'physical' records refer
c to direct I/O while 'logical' records refer to the ACES archive elements.


      external aces_bd_jobarc

c marker(i) : the name of logical record i
c rloc(i)   : the integer index in JOBARC that starts logical record i
c rsize(i)  : the integer-length of logical record i
c nrecs  : the number of physical records in the JOBARC file
c irecwd : the integer-length of a physical record
c irecln : the    recl-length of a physical record

      character*8     marker(1000)
      integer         rloc  (1000),
     &                rsize (1000),
     &                nrecs, irecwd, irecln
      common /jobarc/ marker,
     &                rloc,
     &                rsize,
     &                nrecs, irecwd, irecln
      save   /jobarc/

c bJAUp  : a flag for bombing in get/putrec if aces_ja_init has not been called
c bJAMod : a flag for updating JAINDX in aces_ja_fin

      logical           bJAUp, bJAMod
      common /ja_flags/ bJAUp, bJAMod
      save   /ja_flags/

c jobarc.com : end

c ----------------------------------------------------------------------

c   o make sure the job archive subsystem is up
      if (.not.bJAUp) return

c   o mark the job archive as closed
      call putrec(1,'JOBARC','JADIRTY',1,0)

c   o close the JOBARC file
      inquire(unit=75,opened=bOpened,err=666,iostat=iStat)
      if (bOpened) then
         close(unit=75,status='KEEP',err=666,iostat=iStat)
      end if

c   o only update JAINDX if a record has been added
      if (bJAMod) then

c      o update JAINDX (dump the jobarc common block)
         call gfname('JAINDX',szJAINDX,iJAINDX)
         open(unit=75,file=szJAINDX(1:iJAINDX),
     &        form='UNFORMATTED',status='UNKNOWN',err=666,iostat=iStat)
         rewind(75,err=666,iostat=iStat)
         write(75,err=666,iostat=iStat) marker, rloc, rsize, nrecs
         close(unit=75,status='KEEP',err=666,iostat=iStat)

c      o reset the JOBARC modification flag
         bJAMod = .false.

c     end if (bJAMod)
      end if


c   o turn off the job archive subsystem flag
      bJAUp = .false.

      return

c   o I/O error
 666  print *, '@ACES_JA_FIN: I/O error'
      print '(/)'
      call aces_io_error('ACES_JA_FIN',75,iStat)

c     end subroutine aces_ja_fin
      end

