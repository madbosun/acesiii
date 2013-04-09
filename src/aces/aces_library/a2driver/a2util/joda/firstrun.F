
c This routine performs tasks that should only be executed the first time joda
c is executed. Examples include initializing records that AMEs are expecting to
c exist and possibly checking for files that should not exist.

      subroutine firstrun
      implicit none

      call putrec(1,'JOBARC','JODADONE',1,0)
      call putrec(1,'JOBARC','FNDFDONE',1,1)

      return
      end

