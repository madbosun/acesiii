
cjp for management of multiple moio information for multiple references

      subroutine storemoio
      IMPLICIT INTEGER (A-Z)
#include "bwcc.com"
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &               MOIODS(10,500),MOIOFL(10,500)
      common/io_ptrs/TOTREC(5),TOTWRD(5)

      call icopy(5000, moio  (1,1), 1, bwmoio  (1,1,iref), 1)
      call icopy(5000, moiowd(1,1), 1, bwmoiowd(1,1,iref), 1)
      call icopy(5000, moiosz(1,1), 1, bwmoiosz(1,1,iref), 1)
      call icopy(5000, moiods(1,1), 1, bwmoiods(1,1,iref), 1)
      call icopy(5000, moiofl(1,1), 1, bwmoiofl(1,1,iref), 1)

      call icopy(5,totrec(1),1,bwtotrec(1,iref),1)
      call icopy(5,totwrd(1),1,bwtotwrd(1,iref),1)

      return
      end
