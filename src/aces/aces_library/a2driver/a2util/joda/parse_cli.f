














































































































































































































































































































































































































































































































































      subroutine parse_cli
      implicit none

      integer iLast
      integer i, f_iargc, num_args
      character arg*80
      logical bDoNext

      character*1 achar
      integer  f_strpbrk, c_atol
      external f_strpbrk, c_atol

c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end

c ----------------------------------------------------------------------

      num_args = f_iargc()

c   o loop over arguments
      bDoNext = .true.
      do i = 1, num_args

         if (bDoNext) then
c         o process this argument

            arg=' '
            call f_getarg(i,arg)

c         o MPI_COMM_SIZE
            if (arg(1:6).eq.'-procs') then
               if (i+1.le.num_args) then
                  call f_getarg(i+1,arg)
                  iLast = f_strpbrk(arg,' ')
                  arg(iLast:iLast) = achar(0)
                  nprocs = c_atol(arg)
                  bDoNext = .false.
                  print *, 'nprocs = ',nprocs
               else
                  print *, '@PARSE_CLI: -procs is missing an argument'
                  call aces_exit(1)
               end if
            end if

c         o MPI_COMM_RANK
            if (arg(1:5).eq.'-rank') then
               if (i+1.le.num_args) then
                  call f_getarg(i+1,arg)
                  iLast = f_strpbrk(arg,' ')
                  arg(iLast:iLast) = achar(0)
                  irank = c_atol(arg)
                  bDoNext = .false.
                  print *, 'irank = ',irank
               else
                  print *, '@PARSE_CLI: -rank is missing an argument'
                  call aces_exit(1)
               end if
            end if

c        else if (.not.bDoNext) then
         else

c         o resume processing arguments
            bDoNext = .true.

c        end if (bDoNext)
         end if

c     end do i = 1, num_args
      end do

c ----------------------------------------------------------------------

c   o check consistency

c   o MPI_COMM_SIZE
      if (nprocs.lt.1) then
         print *, '@PARSE_CLI: resetting number of processes to 1'
         nprocs = 1
      end if

c   o MPI_COMM_RANK
      if (irank.lt.0) then
         print *, '@PARSE_CLI: resetting process rank to 0'
         irank = 0
      end if
      if (irank.ge.nprocs) then
         print *, '@PARSE_CLI: resetting process rank to ',nprocs-1
         irank = nprocs-1
      end if

c ----------------------------------------------------------------------

      return
c     end subroutine parse_cli
      end

