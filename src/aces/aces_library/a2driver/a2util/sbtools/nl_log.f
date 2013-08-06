
c This routine reads a logical value.  If the keyword has a value,
c the value must be one of the (case-insensitive) following:
c   yes, true,   1,  on,  y,  t
c   no,  false,  0,  off  n,  f
c If it does not have a value, it is set to .true.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_log(key,def,var)
      implicit none

      character*(*) key
      logical var,def



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



      character*5 val,defval
      logical nl_key,present,streq
      integer strlen

      callstack_prev=callstack_curr
      callstack_curr='NL_LOG'

      var=def

      val=' '
      present=nl_key(key,val)
      if (present) then
         var=.true.
         if (strlen(val).gt.0) then
            if (streq(val,'yes',.true.) .or.
     &          streq(val,'true',.true.) .or.
     &          streq(val,'y',.true.) .or.
     &          streq(val,'t',.true.) .or.
     &          streq(val,'1',.true.) .or.
     &          streq(val,'on',.true.)) then
               var=.true.
            else if (streq(val,'no',.true.) .or.
     &               streq(val,'false',.true.) .or.
     &               streq(val,'n',.true.) .or.
     &               streq(val,'f',.true.) .or.
     &               streq(val,'0',.true.) .or.
     &               streq(val,'off',.true.)) then
               var=.false.
            else
               write(*,*) '@NL_LOG: invalid value for ',key
               call errex
            end if
         end if
      end if

      if (prt_nl .and. .not.streq(key,'print_nl',.true.)) then
         val='false'
         if (var) val=' true'
         defval='false'
         if (def) defval=' true'
         if (val.eq.defval) then
            write(*,900) key,'logical',val
  900       format(a20,2x,a10,2x,a20)
         else
            write(*,910) key,'logical',val,defval
  910       format(a20,2x,a10,2x,a20,2x,a20)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_log
      end

