
cjp for management of multiple moio information for multiple references

      subroutine loadmoio
      IMPLICIT INTEGER (A-Z)
#include "bwcc.com"
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &               MOIODS(10,500),MOIOFL(10,500)
      common/io_ptrs/TOTREC(5),TOTWRD(5)

      call icopy(5000, bwmoio  (1,1,iref), 1, moio  (1,1), 1)
      call icopy(5000, bwmoiowd(1,1,iref), 1, moiowd(1,1), 1)
      call icopy(5000, bwmoiosz(1,1,iref), 1, moiosz(1,1), 1)
      call icopy(5000, bwmoiods(1,1,iref), 1, moiods(1,1), 1)
      call icopy(5000, bwmoiofl(1,1,iref), 1, moiofl(1,1), 1)

      call icopy(5,bwtotrec(1,iref),1,totrec(1),1)
      call icopy(5,bwtotwrd(1,iref),1,totwrd(1),1)

      return
      end
