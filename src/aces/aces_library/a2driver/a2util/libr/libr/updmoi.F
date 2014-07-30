
c This routine was the array (list) creation utility, but that function has been
c transferred to aces_list_touch. This routine is now kept for compatibility.

c NOTE: We cannot get rid of it yet. Calling updmoi with IENTER=1 would attempt
c to reset the free-space pointers of the respective storage file. This action
c is no longer allowed, and although reasonable effort has been made to find any
c residual reset calls, it is possible some may still exist. After some time has
c passed, we can ignore IENTER entirely.

      subroutine updmoi(lstdis,dissiz,lstspn,lstnum,ienter,ioff)
      implicit none

c ARGUMENTS
      integer lstdis, dissiz, lstspn, lstnum, ienter, ioff

c INTERNAL VARIABLES
      integer iTmp

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert ienter=0
      if (ienter.ne.0) then
         print *, '@UPDMOI: Assertion failed.'
         print *, '   ienter = ',ienter
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      call aces_list_touch(lstdis,dissiz,lstspn,lstnum,ioff)

      return
c     end subroutine updmoi
      end

