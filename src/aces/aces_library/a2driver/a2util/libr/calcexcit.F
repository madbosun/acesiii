      subroutine calcexcit(erank,efrom,eto,espin,jr,ir,norb)

      implicit none
#include "bwcc.com"
      integer norb,erank,efrom,eto,espin,jr,ir,i,is,erankfrom,erankto
      dimension efrom(maxexcit), eto(maxexcit),espin(maxexcit)
      dimension norb(2)

      erank=0
      if(ir.eq.jr) return

cjp spin loop
      do is=1,2
      erankfrom=erank
      erankto=erank
cjp find all differences in occupations
      do i=1,norb(is)
      if(iocc(i,ir,is).ne.iocc(i,jr,is)) then
          if(iocc(i,ir,is) .eq. 1 ) then
              erankfrom=erankfrom+1
              efrom(erankfrom) = iphnum(i,ir,is)
          else
              erankto=erankto+1
              eto(erankto) = iphnum(i,ir,is)
          endif
      endif
      enddo

      if(erankfrom.ne.erankto) then
         write(6,*) 'inconsistent number of electrons in refs. ',ir,jr
         stop 'inconsistent number of electrons in references'
      endif
cjp sort indices
      if(erankfrom-erank.gt.1) then
         call sortint(erankfrom-erank,efrom(erank+1))
         call sortint(erankfrom-erank,eto(erank+1))
      endif
cjp store spin
      do i=erank+1,erankfrom
      espin(i)=is
      enddo
cjp increase erank
      erank=erankfrom
cjp
      enddo
cjp spin loop

      return
      end
