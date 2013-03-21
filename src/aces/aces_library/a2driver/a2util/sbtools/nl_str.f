
c This reads a string value.  If it isn't in the namelist,
c it defaults to the value passed in.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_str(key,def,var)
      implicit none

      character*(*) key,def,var



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



      character*(nllinelen) str1,str2
      logical nl_key,present
      integer strlen

      callstack_prev=callstack_curr
      callstack_curr='NL_STR'

      present=nl_key(key,var)
      if (.not.present) var=def

      if (prt_nl) then
         if (strlen(var).le.20) then
            str1=var
         else
            str1=var(1:18)//' $'
         end if
         if (strlen(def).le.20) then
            str2=def
         else
            str2=def(1:18)//' $'
         end if
         if (var(1:strlen(var)).eq.def(1:strlen(def))) then
            write(*,900) key,'string',str1
  900       format(a20,2x,a10,2x,a20)
         else
            write(*,910) key,'string',str1,str2
  910       format(a20,2x,a10,2x,a20,2x,a20)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_str
      end

