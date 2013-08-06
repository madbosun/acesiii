
c This routine prints a summary of all the records in the JOBARC file up to
c the next open slot.

      subroutine aces_ja_summary
      implicit none

c INTERNAL VARIABLES
      integer i
      logical bNotDone

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
         print *, '@ACES_JA_SUMMARY: Assertion failed.'
         print *, '   bJAUp = ',bJAUp
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)

c ----------------------------------------------------------------------

      print '(/)'
      print '(8x,64a)', ('-',i=1,64)
      print '(8x,21x,a22)', 'SUMMARY OF JOB ARCHIVE'
      print '(8x,64a)', ('-',i=1,64)
      print '(8x,a5,8x,3x,a4,3x,8x,a13,8x,a12)',
     &          'INDEX', 'NAME', 'ADDRESS (INT)', 'LENGTH (INT)'
      print '(8x,64a)', ('-',i=1,64)
      i = 1
      bNotDone = .true.
      do while (bNotDone)
      print '(8x,i5,8x,a1,a8,a1,8x,i13,8x,i12)',
     &           i, '"',marker(i),'"', rloc(i), rsize(i)
      if (iand(i,3).eq.0) print '(/)'
      bNotDone = ((marker(i).ne.'OPENSLOT').and.(i.le.1000))
      i = i+1
      end do
      print '(8x,64a)', ('-',i=1,64)
      print '(/)'

      return
c     end subroutine aces_ja_summary
      end

