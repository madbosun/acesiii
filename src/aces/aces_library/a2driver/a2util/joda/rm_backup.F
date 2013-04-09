
c This routine removes the CURRENT backup directory.

#include "flags.h"

      subroutine rm_backup
      implicit none

      character*80 szSaveDir
      integer i, iTmp
#include "jodaflags.com"
      integer ishell
      external ishell

      if (iflags(h_IFLAGS_restart).ne.0) then
         call gfname('SAVEDIR',szSaveDir,i)
         print *, 'Removing ',szSaveDir(1:i),'/CURRENT . . .'
         iTmp=ishell('/bin/rm -rf '//szSaveDir(1:i)//'/CURRENT')
         if (iTmp.eq.0) then
            print *, 'The backup directory was successfully removed.'
         else
            print *, '@RM_BACKUP: The directory ',szSaveDir(1:i),
     &               '/CURRENT could not be removed.'
            print *, '          Future calculations may be ',
     &               'contaminated.'
         end if
c     end if (restart)
      end if

      return
c     end subroutine rm_backup
      end

