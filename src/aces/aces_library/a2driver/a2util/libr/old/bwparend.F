      subroutine bwparend
c finish for BWCC parallelization
      implicit none
#include "bwcc.com"
#ifdef MPI
#include "mpif.h"
       integer ierror
#endif

#ifdef MPI
      write(6,*)' FINAL BARRIER MPI: nproc, myproc ',nproc,myproc
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      write(6,*)' END MPI: nproc, myproc ',nproc,myproc
      call MPI_Finalize(ierror)
#endif
      return
      end
