      subroutine rdcchr(iunit,ncore,igotit)

      implicit double precision (a-h,o-z)
c
c --- <rdcchr> reads the number of core-electrons from
c ---          the ecpdata record
c
c     iunit    : unitnumber of file containing the ecpdata  (input)
c     ncore    : number of core electrons                   (output)
c     igotit   : check of successfull search                (output)
c
      character*80 zeile,scrzl
      character*8 lsymb
      character cdummy
c
      data lsymb /'spdfghij'/

      igotit=-1
c
  100 read(iunit,'(a)',end=900) zeile
        call wisch(zeile,80)
        if(zeile(1:1).eq.' '.or.zeile(1:1).eq.'#') goto 100
c     --- asterisk terminates input
        if(zeile(1:1).eq.'*') goto 900
        icore=index(zeile,'NCORE')
        scrzl=zeile(icore+5:80)
        call rdebbs(scrzl,80,1,0,0,isucc,iwert,rdummy,cdummy)
        if(isucc.le.0) then
          write(6,'(/,2x,a,/)')
     1     'CANNOT READ NCORE = NUMBER OF CORE ELECTRONS'
          return
        else
          igotit=1
        endif
        ncore=iwert
c
      return
c
  900 continue
      write(6,'(/,2x,a,/)') 'RDECP : END OF FILE ENCOUNTERED - SORRY'
      return
c
      end
