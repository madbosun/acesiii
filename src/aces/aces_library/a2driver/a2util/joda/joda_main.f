      subroutine joda_main()

      integer pass1, i
      logical fd, geomopt, num_grad

      integer iMemMin, iMemInc
      parameter (iMemMin=2097152)
      parameter (iMemInc=1048576)

      integer ishell
      external ishell



c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end












c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
      logical ignore
      common /restart_com/ ignore






































































































































































































































































































































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/






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

      geomopt=(iflags2(5).ne.0)
      num_grad=(iflags2(138).eq.1)
      fd = (iflags(54).eq.3.or.(geomopt.and.num_grad))

      if (fd) then

c      o allocate memory for symcor
         icrsiz = iflags(36)
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

         if (ignore) then
c         o mid-stream -> keep going
            call getrec(1,'JOBARC','FNDFDONE',1,i)
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
            else
               call geopt
               if (pass1.eq.0) call putrec(1,'JOBARC','JODADONE',1,1)
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

      return
      end

