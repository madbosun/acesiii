
c initialization routine for BWCC parallelization

      subroutine bwparinit
      implicit none
#include "bwcc.com"
#ifdef MPI
#include "mpif.h"
       integer ierror
#endif

#ifdef MPI
cjp@@@ HANDLING OF STDOUT NOT TO GET MIXED?!!!@@@@@@@
      call mpi_init(ierror)
      call mpi_comm_rank(MPI_COMM_WORLD,myproc,ierror)
      myproc = myproc + 1
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
cYAU      write(*,*) '@BWPARINIT: nproc, myproc ',nproc,myproc
cjp change working directory appropriatelly
cjp first assuming common working directory for all processes
      if (nproc.gt.1) call bwchdir(myproc)
#else
      nproc  = 1
      myproc = 1
#endif

      return
      end
