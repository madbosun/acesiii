
c This routine removes all the records in the JOBARC file at the specified
c record and offset. It does not decrease the size of the JOBARC file, it merely
c resets the lookup table data.

c Example:
c  - erase all records after NEXTGEOM
c    call aces_ja_truncate('NEXTGEOM',1)

c INPUT
c char*(*) SZRECORD : the name of the record to measure the offset from
c int      IOFFSET  : the offset from the record index to truncate at

      subroutine aces_ja_truncate(szRecord,iOffset)
      implicit none

c ARGUMENTS
      character*(*) szRecord
      integer iOffset

c EXTERNAL FUNCTIONS
      integer iszeq

c INTERNAL VARIABLES
      integer ndx, i

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

      i = 0
c   o assert job archive subsystem is up
      if (.not.bJAUp) then
         print *, '@ACES_JA_TRUNCATE: Assertion failed.'
         print *, '   szRecord = ',szRecord
         print *, '   bJAUp    = ',bJAUp
         i = 1
      end if
c   o the record name must be between 1 and 8 characters
      if ((len(szRecord).lt.1).or.(len(szRecord).gt.8)) then
         print *, '@ACES_JA_TRUNCATE: Assertion failed.'
         print *, '   szRecord = ',szRecord
         i = 1
      end if
c   o the record cannot be named 'OPENSLOT'
      if (len(szRecord).eq.8) then
      if (szRecord.eq.'OPENSLOT') then
         print *, '@ACES_JA_TRUNCATE: Assertion failed.'
         print *, '   szRecord = ',szRecord
         i = 1
      end if
      end if
      if (i.ne.0) call aces_exit(i)

c ----------------------------------------------------------------------

c   o get the index of the requested record
      ndx = iszeq(1000,marker,1,szRecord)

c   o the record is not found or is at the last position
      if ((ndx.eq.0).or.(ndx.eq.1000)) return

c   o point to the truncation record
      ndx = ndx + iOffset

c   o the position is out of bounds
      if ((ndx.lt.1).or.(ndx.gt.1000)) return

c   o the record is already empty
      if (marker(ndx).eq.'OPENSLOT') return

c   o reset the remaining structures
      do i = ndx, 1000
         marker(i) = 'OPENSLOT'
      end do
      call izero(rloc(ndx), 1+1000-ndx) 
      call izero(rsize(ndx),1+1000-ndx) 

c   o mark JOBARC as modified
      bJAMod = .true.

      return
c     end subroutine aces_ja_truncate
      end

