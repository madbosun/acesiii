
c This routine defines the default values for the common blocks in
c parallel_aces.com.

      blockdata aces_bd_parallel_aces
      implicit none
c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.





      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end
      data nprocs,irank /1,0/
      data szcpuname /'localhost'/
      data  icpuname /9/
      end

