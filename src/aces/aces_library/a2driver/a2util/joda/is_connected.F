
      function is_connected(atom1, atom2, connectiontable, natoms,
     &                      connectlist)
C
C Determines if two atoms are connected by constructing a full
C connection list and then searching it.
C
      implicit none
C
C Input Variables
C
      integer atom1,atom2,natoms,connectiontable(natoms,natoms)
C
C Function Variables
C
      logical is_connected
C
C Pre-allocated Local Variables
C
      integer connectlist(natoms), i
C
C External Procedures
C
      external followconnect

C Local Variables

      integer listatoms
      logical shortcircuit
C - - -- - - - - - - - - - - - - -- - - - - - - -- - - -  -- - 

      if (connectiontable(atom1,atom2).eq.1) then
         is_connected=.true.
         return
      else
         call izero(connectlist,natoms)
C
C Follow atom1 connections
C
         listatoms=0
         shortcircuit=.false.
         call followconnect(atom2,atom1,connectlist,listatoms,
     &      connectiontable,natoms,followconnect,shortcircuit)
         is_connected=shortcircuit
      endif
C
#ifdef _DEBUG_LVLM1
      Print*, "The connectionlist in is_connected"
      Print*, "atom1 & atom2", (connectlist(i), i=1, natoms),
     &         atom1, atom2
#endif
c
      return
      end
C
      subroutine followconnect(atom,taratom,list,listatoms,
     &   connectiontable,natoms,follow2,shortcircuit)
C
C Constructs list of connections pseudo-recursively
C by calling follow2 (which is a pointer to followconnect)
C
      implicit none
C Input Variables
      integer atom,natoms,connectiontable(natoms,natoms),taratom

C Input/Output Variables
      integer list(natoms),listatoms
      logical shortcircuit

C External Procedures
      external follow2,bin_search
      logical bin_search

C Local Variables
      integer ii
      integer iatom 
      logical new
C - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - - -- - -
      
      if (atom.eq.taratom) then
         shortcircuit=.true.
      else
C
C Determine if atom is already in the list
C
         if (listatoms.gt.0) then
            new=.not.bin_search(atom,list,listatoms)
         else
            new=.true.
         endif
      
         if (new) then
C Insert in sorted list
            ii=listatoms
            do while ((ii.ge.1).and.(list(ii).gt.atom))
               list(ii+1)=list(ii)
               ii=ii-1
            end do
            list(ii+1)=atom
            listatoms=listatoms+1
C Follow the connections down a level
            do 2 ii=1,natoms
               if (connectiontable(ii,atom).eq.1) then
                  iatom = ii
                  call follow2(iatom,taratom,list,listatoms,
     &               connectiontable,
     &               natoms,follow2,shortcircuit)
C Short circuit indicates that the target atom taratom has been found
                  if (shortcircuit) return
               endif
 2          continue
         endif
      endif
      
      return
      end

      function bin_search(atom,list,listatoms)
C
C Binary search routine returns true if atom is in list false 
C otherwise
C
      implicit none
C Input Variables
      integer atom,listatoms,list(listatoms)

C Function Variables
      logical bin_search

C Local Variables
      integer left,right,mid
C - - - - -- - -- - - - - - -- - - - - -- - - - - - - - - - - -
      left=1
      right=listatoms
      do while (left .le. right)
         mid=int(dble(right-left)/2.0d0)+left
         if (list(mid).eq.atom) then
            bin_search=.true.
            return
         else if (atom .lt. list(mid)) then
            right=mid-1
         else
            left=mid+1
         endif
      end do
      bin_search=.false.
      return
      end
