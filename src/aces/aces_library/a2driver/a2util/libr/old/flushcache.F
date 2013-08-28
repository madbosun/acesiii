      SUBROUTINE flushcache
C
cjp
cjp   necessary for work with iref-replicated moints, moabcd files
cjp
      IMPLICIT INTEGER (A-Z)
      DIMENSION ZLIST(1)
      COMMON // ICORE(1)
      COMMON /IOPOS/ ICRSIZ,ICHCSZ,IOFF(2),LENREC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /MACHSP2/ MASK1,MASK2,ISHFSZ
      COMMON /CACHEINF/ CACHNUM,CACHNMP1,CACHDIR(100),CACHPOS(100),
     &                  CACHFILE(100),CACHMOD(100),OLDEST
      COMMON /FILSPC/ ILNBUF,IPRCLN,IPRCWD
C
cjp for debug purposes in order to identify exactly, when some
cjp modifications of files occured, switch off caching mechanism
cjp by flushing caches immediatelly

cjp will be needed anyway in the mr-cc to flush after work with a set of amplitudes
cjp of one reference before starting with the other one

      IPACK(I,J)=IOR(J,ISHFT(I-49,ISHFSZ))
      UPACKR(I) =IAND(I,MASK1)
      UPACKF(I) =IAND(ISHFT(I,-ISHFSZ),MASK2)+49

cjp this flush procedure is taken from termio
      do 15 i=1,cachnum
       file=cachfile(i)
       if(file.ne.0)then
        record=upackr(cachdir(i))
        ipos  =cachpos(i)
        if(cachmod(i).ne.0.and.record.ne.0)
     +       call aces_io_write(file,record,icore(ipos),lenrec)
        cachmod(i)=0
       endif
15    continue

cjp mark also the cache positions as empty
cjp in order that subsequent calls of getlst/fetch do not hit the cache
cjp full of old stuff from other reference configuration!!!
cjp code from initio
      call izero(cachdir,100)
      call izero(cachfile,100)
      call izero(cachmod,100)
      oldest=1
      return
      end
