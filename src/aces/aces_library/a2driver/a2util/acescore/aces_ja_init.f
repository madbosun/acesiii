
c This routine initializes the job archive environment (via the JAINDX file)
c in order for getrec and putrec to work properly.







      subroutine aces_ja_init
      implicit none

c INTERNAL VARIABLES
      character*80 szJOBARC, szJAINDX
      integer       iJOBARC,  iJAINDX
      integer iBuf(128), i, iStat
      logical bExist, bOpened

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
c icdacc.com : begin
c Nevin 8-30-95 added record length common to facilitate change from
c Bytes to Words for SGI and DecAlpha
      integer         idaccm
      common /icdacc/ idaccm
c icdacc.com : end

c ----------------------------------------------------------------------

c   o die instead of silently return if the job archive subsystem is already up
      if (bJAUp) then
         print *, '@ACES_JA_INIT: The job archive subsystem is already',
     &            ' initialized.'
         call aces_exit(1)
      end if

c   o turn on the job archive subsystem flag
      bJAUp = .true.

c   o get the external JOBARC file name
      call gfname('JOBARC',szJOBARC,iJOBARC)

c   o initialize the jobarc common block
      call gfname('JAINDX',szJAINDX,iJAINDX)
      inquire(file=szJAINDX(1:iJAINDX),exist=bExist,
     &        err=666,iostat=iStat)
      if (bExist) then
c      o JOBARC had better exist
         inquire(file=szJOBARC(1:iJOBARC),exist=bExist,
     &           err=666,iostat=iStat)
         if (.not.bExist) then
            print *, '@ACES_JA_INIT: JOBARC does not exist'
            call aces_exit(1)
         end if
c      o process JAINDX
         open(unit=75,file=szJAINDX(1:iJAINDX),
     &        form='UNFORMATTED',status='OLD',err=666,iostat=iStat)
         rewind(75,err=666,iostat=iStat)
         read(75,err=666,iostat=iStat) marker, rloc, rsize, nrecs
         close(unit=75,status='KEEP',err=666,iostat=iStat)
      else
c      o reset the records
         do i = 1, 1000
            marker(i) = 'OPENSLOT'
         end do
         rloc(1) = 1
         nrecs   = 0
      end if
      irecwd = 128
      irecln = 128*idaccm

c   o condition the JOBARC file
      inquire(file=szJOBARC(1:iJOBARC),exist=bExist,opened=bOpened,
     &        err=666,iostat=iStat)
      if (bExist) then
c      o open the JOBARC file
         if (.not.bOpened) then
            open(unit=75,file=szJOBARC(1:iJOBARC),
     &           form='UNFORMATTED',access='DIRECT',
     &           status='OLD',recl=irecln,err=666,iostat=iStat)
         end if
      else
c      o write out one all-zero record (bug fix: Ajith 03/25/97)
         call izero(iBuf,128)
         open(unit=75,file=szJOBARC(1:iJOBARC),
     &        form='UNFORMATTED',access='DIRECT',
     &        status='NEW',recl=irecln,err=666,iostat=iStat)
         write(unit=75,rec=1,err=666,iostat=iStat) iBuf
      end if

c   o check if the job archive was properly closed and mark it as open
      call getrec(-1,'JOBARC','JADIRTY',1,i)
      if (i.ne.0) then
         print '(/)'
         print *, '@ACES_JA_INIT: WARNING - The job archive was not ',
     &            'finalized by the previous'
         print *, '               ACES Member Executable. Any records ',
     &            'added by that process have'
         print *, '               been lost.'
         print '(/)'
      else
         call putrec(1,'JOBARC','JADIRTY',1,1)
      end if

      return

c   o I/O error
 666  print *, '@ACES_JA_INIT: I/O error'
      print '(/)'
      call aces_io_error('ACES_JA_INIT',75,iStat)

c     end subroutine aces_ja_init
      end

