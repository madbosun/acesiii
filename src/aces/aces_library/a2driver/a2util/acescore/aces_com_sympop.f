
c This routine initializes the sympop common block, but isytyp is loaded in
c aces_io_init.

      subroutine aces_com_sympop
      implicit none

c INTERNAL VARIABLES
      integer i
      character*8 symlst(22)
      data symlst /'SVAVA0X', 'SVBVB0X',
     &             'SOAOA0X', 'SOBOB0X',
     &             'SVAVA1X', 'SVBVB1X',
     &             'SOAOA1X', 'SOBOB1X',
     &             'SVAOA2X', 'SVBOB2X',
     &             'SOBVA2X', 'SVBOA2X',
     &             'SVAVB2X', 'SOAOB2X',
     &             'SVAVB2X', 'SOAVA2X',
     &             'SOBVB2X', 'SOAVB2X',
     &             'SVAVA2X', 'SVBVB2X',
     &             'SOAOA2X', 'SOBOB2X'/

c COMMON BLOCKS
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

c ----------------------------------------------------------------------

      do i = 1, 22
         call getrec(-1,'JOBARC',symlst(i),nirrep,irpdpd(1,i))
      end do

      return
c     end subroutine aces_com_sympop
      end

