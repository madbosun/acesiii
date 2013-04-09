      subroutine rdebbs(scr,n,kint,kreal,kchar,isucc,iwert,rwert,cwert)
C
      implicit double precision (a-h,o-z)
c
c     --- this routine reads one integer/real/character variable
c         as specified by the mask kint/kreal/kchar (two values of
c         which have to be zero) from the input string <scr>
c         it is assumed that integer/real numbers do not need
c         more than 32 characters of space (if wanted, pump that up)
c         if the mask value is positive, the input characters
c         will be cleaned after reading if the mask value is
c         positive
c
c         NOTE : abs(kchar) = length of the output string <cwert>
c
      integer kint,kreal,kchar
      integer isucc,iwert
      double precision rwert
      character*(*) cwert,scr
      character*32 readit
      logical prtout
c
      prtout=n.lt.0
c
      if(prtout) write(6,601) scr
  601 format(/,' input string = ',/,a,/)
c
c     --- if read succeeds, isucc contains the length of the string
c
c     --- check input parameters
      if((kint.ne.0.and.kreal.ne.0).or.
     1   (kint.ne.0.and.kchar.ne.0).or.
     2   (kreal.ne.0.and.kchar.ne.0)) stop ' abuse of rdebbs '
c
      iwert=0
      rwert=0.d0
      cwert=' '
c
      isucc=iblank(scr)-1
c
      if(isucc.eq.0) return
c
      if((kint.ne.0.or.kreal.ne.0).and.isucc.gt.32) then
        write(6,901)
  901   format(/,' integer/real number with more than 32 digits ',/)
        goto 800
      elseif(kchar.ne.0.and.abs(kchar).lt.isucc) then
        write(6,902) iabs(kchar)
  902       format(/,' i/o-error : input string is longer than ',i3,
     1             ' characters ',/)
        goto 800
      endif
c
      nbl=32-isucc
      if(kint.ne.0) then
        iact=kint
        if(nbl.gt.0) readit(1:nbl)=' '
        readit(nbl+1:32)=scr(1:isucc)
        read(readit,'(i32)',err=701) iwert
      elseif(kreal.ne.0) then
        iact=kreal
        if(nbl.gt.0) readit(1:nbl)=' '
        readit(nbl+1:32)=scr(1:isucc)
        read(readit,'(g32.0)',err=702) rwert
      elseif(kchar.ne.0) then
        iact=kchar
        cwert=scr(1:isucc) 
      endif
c     --- delete the string if iact>0 (iact<0 allows second read)
      if(iact.gt.0) scr(1:isucc)=' '
c
      return
c
  701 write(6,751)
  751 format(/,' i/o-error : input variable is not integer ',/)
      goto 800
  702 write(6,752)
  752 format(/,' i/o-error : input variable is not real ',/)
c
  800 write(6,991) scr
  991 format(/,
     1 ' WARNING : <rdebbs> could not read properly from string ',/,a,/)
      isucc=-1
c
      return
      end
