
c This routine opens the ZMAT file and searches for a namelist.
c If found, it is read into the array namelist for parsing.
c
c If printdef is true, the namelist will default to printing (though it
c can be overridden by the PRINT_NL keyword in the namelist itself).

      subroutine nl_init(title,err,printdef)
      implicit none

      character*(*) title
      integer err
      logical printdef



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



      character*(nllinelen) line, szTitle
      integer zio, linenum, l, strlen
      logical streq, bPrint, bOpen

      callstack_prev = callstack_curr
      callstack_curr = 'NL_INIT'

      szTitle = '*'//title
      err     = 0
      bPrint  = printdef

c Open ZMAT file and search for the namelist.
      zio = 10
      inquire(unit=zio,opened=bOpen)
      do while (bOpen)
         zio = zio+1
         inquire(unit=zio,opened=bOpen)
      end do
      open(unit=zio,file='ZMAT',form='formatted',err=999,iostat=err)

      linenum=0
      nlnumline=0
      l=strlen(szTitle)
   10 continue
         read(zio,'(a)',end=999) line
         linenum=linenum+1
      if (strlen(line).eq.0 .or.
     &    line(1:1).eq.'#' .or.
     &    .not.streq(szTitle,line(1:l),.true.)
     &   ) goto 10

c Once the namelist is found, remove the '*TITLE' from it and if
c anything is left, put it as the first line of the namelist text
      line=line(l+1:strlen(line))
      if (strlen(line).gt.0) then
         call upcase(line)
         nltext(1)=line
         nlnumline=1
      end if

c Read until we are done with the file or find a line starting with '*'
      line=' '
      do while (line(1:1).ne.'*')
         read(zio,'(a)',end=999) line
         linenum=linenum+1
         if (strlen(line).ne.0 .and.
     &       line(1:1).ne.'#' .and.
     &       line(1:1).ne.'*') then
            if (nlnumline.eq.nlmaxline) then
               write(*,'(a)')
     &            '@NL_INIT: maximum length of namelist exceeded'
               call errex
            end if
            nlnumline=nlnumline+1
            call upcase(line)
            nltext(nlnumline)=line
         end if
      end do

  999 continue
      close(zio)
      call nl_log('print_nl',bPrint,prt_nl)
      if (prt_nl) call nl_prttop(title)
      if (nlnumline.eq.0) err=1

      callstack_curr = callstack_prev
      return
c     end subroutine nl_init
      end

