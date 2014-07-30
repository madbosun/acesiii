
c This routine checks to see if any invalid keywords are left in line.

      subroutine nl_term
      implicit none


c This common block contains information about a namelist which is being
c parsed.

c   nlmaxline   : the maximum number of lines which can be in the namelist
c   nllinelen   : the maximum length of each line

c   nlnumline   : the number of lines in the namelist (blank lines are
c                 removed)
c   nltext      : the text in the namelist
c   prt_nl      : a logical which is read in from the namelist to see
c                 if the values read in are printed or not
c
      integer nlmaxline, nllinelen
      parameter(nlmaxline=64)
      parameter(nllinelen=132)

      character*(nllinelen) nltext(nlmaxline)
      integer nlnumline
      logical prt_nl

      common /namelistc/ nltext
      common /namelist/  nlnumline
      common /namelistl/ prt_nl
      save /namelist/
      save /namelistc/
      save /namelistl/




c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/


      integer i
      callstack_prev=callstack_curr
      callstack_curr='NL_TERM'
      if (prt_nl) call nl_prtbot
      if (nlnumline.ne.0) then
         write(*,*) '@NL_TERM: invalid namelist entries'
         write(*,*) '  Errors somewhere around '
         write(*,*)
         do i = 1, nlnumline
            write(*,*) nltext(i)
         end do
         call errex
      end if
      callstack_curr=callstack_prev
      return
      end

