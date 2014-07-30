
c This routine reads an integer value.  If it is not in the namelist,
c it defaults to the value passed in.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_int(key,def,var)
      implicit none

      character*(*) key
      integer def,var



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



      character*(nllinelen) val
      logical nl_key,present
      integer err

      callstack_prev=callstack_curr
      callstack_curr='NL_INT'

      var=def

      val=' '
      present=nl_key(key,val)
      if (present) then
         call str2int(val,var,err)
         if (err.ne.0) then
            write(*,*) '@NL_INT: invalid integer for ',key
            call errex
         end if
      end if

      if (prt_nl) then
         if (var.eq.def) then
            write(*,900) key,'integer',var
  900       format(a20,2x,a10,2x,i20)
         else
            write(*,910) key,'integer',var,def
  910       format(a20,2x,a10,2x,i20,2x,i20)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_int
      end

