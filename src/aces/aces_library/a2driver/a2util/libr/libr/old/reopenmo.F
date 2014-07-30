      subroutine reopenmo
cjp closes moints and moabcd
cjp abd reopens them with new iref
cjp
      implicit integer (a-z)
      LOGICAL YESNO,GAMLAM
      CHARACTER*8 NAMES(5)
      CHARACTER*80 FNAME
      INTEGER USEC,USECIN,TSTART,WSTART
      COMMON / / ICORE(1)
      COMMON /INCORE/ ICREC(2),USEC,IXT(2),IMOD(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIODS(10,500),
     &               MOIOSZ(10,500),MOIOFL(10,500)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /IOPOS/ ICRSIZ,ICHCSZ,IOFF(2),LENREC
      COMMON /FILSPC/ IBUFLN,IPRCLN,IPRCWD
      COMMON /CACHEINF/ CACHSIZ,CACHSZP1,CACHDIR(100),CACHPOS(100),
     &                  CACHFILE(100),CACHMOD(100),OLDEST
      DATA NAMES /'MOINTS  ','GAMLAM  ','MOABCD  ','DERINT  ',
     &            'DERGAM  '/
#include "bwcc.com"
c
      if(nref.gt.1) call loadmoio
      do i= 50,52,2
      close(i)
      call gfname(names(i-49),fname,ilength)
      call aces_io_open(i-49)
      enddo
      return
      end
