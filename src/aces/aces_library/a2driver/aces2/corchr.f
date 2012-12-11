      subroutine corchr(cnick,iatnr,nelecp,natoms)

      implicit double precision (a-h,o-z)
c
c --- <corchr> reads the number of core-electrons from the
c ---          file ECPDATA
c
c     cnick    : nicknames of ecp's        (input)
c     iatr     : atomic number of atoms    (input)
c     nelecp   : number of core electrons  (output)
c     natoms   : number of atoms in ZMAT   (input)
c
      character*80 filnam,zeile
      character*80 cnick(natoms)
      dimension nelecp(natoms),iatnr(natoms)
c
c-----------------------------------------------------------------------
c     ecp data:
c-----------------------------------------------------------------------
      iunit=99
      filnam='ECPDATA'
c-----------------------------------------------------------------------
c
c     --- now we are ready for reading
c
      lenfil=index(filnam,' ')-1
      open(unit=iunit,file=filnam(1:lenfil),form='formatted')
      do 400 iecp=1,natoms
        call getstr(iunit,1,igotit,zeile,cnick(iecp))
c --- check if atom has no ecp!
        indi=index(cnick(iecp),'NONE')
c --- or is a dummyatom!
        if (indi.ge.1.or.iatnr(iecp).lt.1) then
          nelecp(iecp)=0
          goto 400
        endif              
c
        if(igotit.le.0) then
          write(6,4091) cnick(iecp),filnam(1:lenfil)
 4091     format(/,1x,' ECPGET ABEND : ECP SET WITH NICKNAME ',A,/,
     1             1x,' IS NOT CONTAINED WITHIN FILE ',A,/)
          goto 999
        endif
c
  11    call getstr(iunit,-1,ifound,zeile,'*')
c     '*' has to be the first character in the line :
        if(ifound.gt.1) goto 11
        if(ifound.gt.0) call rdcchr(iunit,nelecp(iecp),igotit)
c
        if(igotit.le.0.or.ifound.le.0) then
          write(6,4092) cnick(iecp)
 4092     format(/,1x,' ecpget abend : ecp set with nickname ',a,/,
     1             1x,' cannot be read in properly ',/)
          goto 999
        endif
c
  400 continue
c
      close(iunit)
      return
c
  999 continue
      write(6,901)
  901 format(/,1x,' <CORCHR> : CANNOT ASSIGN ECP DATA ',/)
      close(iunit)
      CALL ERREX
      end
