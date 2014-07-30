
C Read the Z-matrix and find errors.




      subroutine debug_zmat
      implicit none

c "parameters"
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
c     Maximum string length of terminal lines
      INTEGER LINELEN
      PARAMETER (LINELEN=80)

      character*(linelen) zline, dline
      integer izl(2,7), ref(3)
      logical not_done

      integer atoi, linblnk, i, j, k, iatom

      character*1 achar, czTab, czSpace, czPercent, czComment

c ----------------------------------------------------------------------

      czTab     = achar(9)
      czSpace   = achar(32)
      czComment = achar(35)
      czPercent = achar(37)

      inquire(file=zfil,opened=not_done)
      if (.not.not_done)
     &   open(unit=luz,file=zfil,form='formatted',status='old')
      rewind luz

c   o skip the header (zline=TITLE on end do)
      not_done = .true.
      do while (not_done)
         read(luz,'(a)') zline
         call parsez(zline,izl)
         i = izl(1,1)
         not_done = (i.eq.0).or.(zline(i:i).eq.czPercent)
      end do

c   o first element
      read(luz,'(a)') zline
      call parsez(zline,izl)
      if (izl(1,1).eq.0) then
         write(*,*)
     &      '@DEBUG_ZMAT: There are no elements in the Z-matrix.'
         call errex
      end if
      if (izl(1,2).ne.0) then
         write(*,*)
     &      '@DEBUG_ZMAT: There are extra tokens on the first line.'
         write(*,*)
     &      '             (Perhaps the title line is missing.)'
         write(*,*) ' 1 : ',zline(izl(1,1):izl(2,2)),' ...'
         call errex
      end if

c   o second element
      read(luz,'(a)') zline
      call parsez(zline,izl)
      if (izl(1,3).eq.0) then
         write(*,*)
     &      '@DEBUG_ZMAT: Line 2 must have three values.'
         i = max( izl(2,2) , izl(2,1) )
         write(*,*) ' 2 : ',zline(izl(1,1):i)
         call errex
      end if
      if (izl(1,4).ne.0) then
         write(*,*)
     &      '@DEBUG_ZMAT: There are extra tokens on the second line.'
         write(*,*) ' 2 : ',zline(izl(1,1):izl(2,4)),' ...'
         call errex
      end if

c   o third element
      read(luz,'(a)') zline
      call parsez(zline,izl)
      if (izl(1,5).eq.0) then
         write(*,*)
     &      '@DEBUG_ZMAT: Line 3 must have five values.'
         i = max( izl(2,4) , izl(2,3) , izl(2,2) , izl(2,1) )
         write(*,*) ' 3 : ',zline(izl(1,1):i)
         call errex
      end if
      if (izl(1,6).ne.0) then
         write(*,*)
     &      '@DEBUG_ZMAT: There are extra tokens on the third line.'
         write(*,*) ' 3 : ',zline(izl(1,1):izl(2,6)),' ...'
         call errex
      end if

c   o remaining elements
      iatom = 4
      read(luz,'(a)') zline
      call parsez(zline,izl)
      not_done = .true.
      do while (not_done)

         if (izl(1,7).eq.0) then
            write(*,*)
     &         '@DEBUG_ZMAT: Line ',iatom,' must have seven values.'
            i = max( izl(2,6) , izl(2,5) , izl(2,4) ,
     &               izl(2,3) , izl(2,2) , izl(2,1) )
            write(*,*) iatom,' : ',zline(izl(1,1):i)
            call errex
         end if

         ref(1) = atoi( zline(izl(1,2):izl(2,2)) )
         ref(2) = atoi( zline(izl(1,4):izl(2,4)) )
         ref(3) = atoi( zline(izl(1,6):izl(2,6)) )

c        REUSED REFS

         if (ref(1).eq.ref(2)) then
            dline = zline
            call prep_dline_2(dline,izl,2,4)
            write(*,*) '@DEBUG_ZMAT: References cannot be reused.'
            write(*,*) iatom,' : ',zline(izl(1,1):izl(2,7))
            write(*,*) iatom,' : ',dline(izl(1,1):izl(2,7))
            write(*,*)
         end if

         if (ref(1).eq.ref(3)) then
            dline = zline
            call prep_dline_2(dline,izl,2,6)
            write(*,*) '@DEBUG_ZMAT: References cannot be reused.'
            write(*,*) iatom,' : ',zline(izl(1,1):izl(2,7))
            write(*,*) iatom,' : ',dline(izl(1,1):izl(2,7))
            write(*,*)
         end if

         if (ref(2).eq.ref(3)) then
            dline = zline
            call prep_dline_2(dline,izl,4,6)
            write(*,*) '@DEBUG_ZMAT: References cannot be reused.'
            write(*,*) iatom,' : ',zline(izl(1,1):izl(2,7))
            write(*,*) iatom,' : ',dline(izl(1,1):izl(2,7))
            write(*,*)
         end if

c        ZERO REFS

         do j = 1, 3
            if (ref(j).eq.0) then
               dline = zline
               i = ishft(j,1)
               call prep_dline_2(dline,izl,i,i)
               write(*,*) '@DEBUG_ZMAT: Reference ',j,' for element ',
     &                    iatom,' does not exist.'
               write(*,*) iatom,' : ',zline(izl(1,1):izl(2,7))
               write(*,*) iatom,' : ',dline(izl(1,1):izl(2,7))
               write(*,*)
            end if
         end do

c        UNDEFINED REFS

         do j = 1, 3
            if (ref(j).gt.iatom) then
               dline = zline
               i = ishft(j,1)
               call prep_dline_2(dline,izl,i,i)
               write(*,*) '@DEBUG_ZMAT: Reference ',j,' for element ',
     &                    iatom,' has not yet been defined.'
               write(*,*) iatom,' : ',zline(izl(1,1):izl(2,7))
               write(*,*) iatom,' : ',dline(izl(1,1):izl(2,7))
               write(*,*)
            end if
         end do

c        SELF-REF

         do j = 1, 3
            if (ref(j).eq.iatom) then
               dline = zline
               i = ishft(j,1)
               call prep_dline_2(dline,izl,i,i)
               write(*,*) '@DEBUG_ZMAT: Element ',iatom,
     &                    ' references itself.'
               write(*,*) iatom,' : ',zline(izl(1,1):izl(2,7))
               write(*,*) iatom,' : ',dline(izl(1,1):izl(2,7))
               write(*,*)
            end if
         end do

         iatom = iatom + 1
         read(luz,'(a)') zline
         call parsez(zline,izl)
         not_done = (linblnk(zline(1:80)).ne.0)

c     end do while (not_done)
      end do

      call errex
      return
      end

