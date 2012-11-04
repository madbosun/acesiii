
#include "flags.h"

      program joda
      implicit none

#include "fnamelen.par" 
      character*(fnamelen)fname
      integer pass1, i, i_havegeom, ipost_vib, ilength
      logical fd, geomopt, num_grad, vib_frqs, we_havegeom,
     &        optarc_presnt

      integer iMemMin, iMemInc, Popt_num_fd, do_pes_scan
      parameter (iMemMin=2097152)
      parameter (iMemInc=1048576)

      integer ishell
      external ishell

#include "icore.com"
#include "istart.com"
      logical ignore
      common /restart_com/ ignore
#include "jodaflags.com"

#ifdef _TRAPFPE
c   o since JODA does not call crapsi, have it turn on FPEs
      call trapfpe
#endif /* _TRAPFPE */

c   o initialize the runtime environment
      call aces_init_rte

c   o gather parallel statistics (needed by gfname in ja_init)
      call aces_com_parallel_aces

c   o parse command line (overwrites rank and number of processes)
      call parse_cli

c   o evaluate and repair the health of the current file set
      call dfiles(ignore)

c   o initialize the job archive subsystem
      call aces_ja_init
c
CSSS      call init_coord_com 

      if (.not.ignore) then
c      o this is the first joda run
         call putrec(1,'JOBARC','FIRSTRUN',1,1)
         call putrec(1,'JOBARC','DIRTYFLG',1,0)
         call putrec(1,'JOBARC','JODADONE',1,0)
         call putrec(1,'JOBARC','FNDFDONE',1,1)
         call geopt
      else
c      o load flags
         call getrec(1,'JOBARC','IFLAGS', 100,iflags)
         call getrec(1,'JOBARC','IFLAGS2',500,iflags2)
      end if
c
      call getrec(-1,'JOBARC', 'HAVEGEOM', 1, i_havegeom)
      call getrec(-1,'JOBARC', 'POSTOPFD', 1, popt_num_fd)
      call getrec(-1,'JOBARC', 'PES_SCAN', 1, do_pes_scan)
c
      call gfname('OPTARC', fname, ilength)
      inquire(file=fname(1:ilength), exist=optarc_presnt)
c
      If (i_havegeom .eq. 1) we_havegeom = .true.
      geomopt=(iflags2(h_IFLAGS2_geom_opt).ne.0) .and.
     &         .not. we_havegeom
      num_grad=(iflags2(h_IFLAGS2_grad_calc).eq.2)
      vib_frqs=(iflags(h_IFLAGS_vib).eq.3)
c
      fd = ((iflags(h_IFLAGS_vib).eq.3 .and. we_havegeom)
     &     .or. num_grad .or. (geomopt.and.num_grad))
      if (do_pes_scan .eq. 1 .and. fd) ignore = .false.
c
#ifdef _DEBUG_LVL0
      Print*, "Finite diffs; ignore and geomopt vars:",
     &         ignore, geomopt
      Print*, "The vib calc. related varrs:"
      Print*, "iflags(h_IFLAGS_vib):", iflags(h_IFLAGS_vib)
      Print*, "The finite diffs:", fd
      Print*, "iflags2(h_IFLAGS2_geom_opt):",
     &         iflags2(h_IFLAGS2_geom_opt)
      Print*, "first run of popt num. frq:", popt_num_fd
 
#endif
c
      if (fd) then

c      o allocate memory for symcor
         icrsiz = iflags(h_IFLAGS_mem)
         iCore(1) = 0
         do while ((iCore(1).eq.0).and.(icrsiz.gt.iMemMin))
            call aces_malloc(icrsiz,iCore,i0)
            if (iCore(1).eq.0) icrsiz = icrsiz - iMemInc
         end do
         if (iCore(1).eq.0) then
            print *, '@JODA: unable to allocate at least ',
     &               iMemMin,' integers of memory'
            call aces_exit(1)
         end if
c
          if (ignore) then 
c         o mid-stream -> keep going
            call getrec(1,'JOBARC','FNDFDONE',1,i)
            if (popt_num_fd .eq. 1) call post_opt_update
            if (i.eq.0) call symcor(icore(i0),icrsiz)
         else
c         o first run -> reset finite difference series
            if (geomopt) i=ishell('cp OPTARC OPTARCBK')
            call putrec(1,'JOBARC','FNDFDONE',1,0)
            call symcor(icore(i0),icrsiz)
            ignore=.true.
         end if

         call getrec(1,'JOBARC','PASS1',1,pass1)
         if (pass1.ne.-1) then
c         o vib freqs w/ an grads -OR- geom opts w/ num grads
            call putrec(1,'JOBARC','FIRSTRUN',1,0)
            call putrec(1,'JOBARC','DIRTYFLG',1,1)
            if (geomopt.and.pass1.eq.0) then
               i=ishell('cp OPTARCBK OPTARC')
               call geopt
               call getrec(1,'JOBARC','JODADONE',1,i)
               if (i.ne.1) then
c               o new geom -> reset finite difference series
                  i=ishell('cp OPTARC OPTARCBK')
                  call putrec(1,'JOBARC','FNDFDONE',1,0)
                  call symcor(icore(i0),icrsiz)
                  call geopt
               end if
            else if (vib_frqs .or. num_grad) then 
               if (we_havegeom .and. optarc_presnt) 
     &             i=ishell('rm OPTARC')
               call getrec(1,'JOBARC','FNDFDONE',1,i)
               call geopt
               if (pass1.eq.0) call putrec(1,'JOBARC','JODADONE',1,1)
               call getrec(1,'JOBARC','FNDFDONE',1,i) 
            end if
         end if

      else

         if (ignore) then
c         o joda has run before
            call putrec(1,'JOBARC','FIRSTRUN',1,0)
            call putrec(1,'JOBARC','DIRTYFLG',1,1)
            call geopt
          end if

      end if

c   o finalize the job archive subsystem
      call aces_ja_fin

      end

