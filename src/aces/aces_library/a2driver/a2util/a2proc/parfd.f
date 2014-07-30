
c This routine prepares or loads data needed for parallel finite differences.
c The goto routine for what should be done is symcor/upd_fd.F












































































































































































































































































































































































































































































































































      subroutine parfd(args,dimargs)
      implicit none

c ARGUMENT LIST
      integer dimargs
      character*80 args(*)

c INTERNAL VARIABLES
      double precision dDipXYZ(3), dPolXYZ(3,3), dTmp
      integer nAtom, nSize, nPoint, iOff, iLast, iTmp
      integer lFree, lGeom, lGrdXYZ, lGrdInt, lIMAP
      character*80 szRecName
      character*4 szType
      logical bUpdate, bDump, bAnGrad, bExist, bNotDone

c EXTERNAL FUNCTIONS
      integer  linblnk, atoi
      external linblnk, atoi

c COMMON BLOCKS


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end





c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end

c ----------------------------------------------------------------------

c   o load system size and number of displacements
      call getrec(1,'JOBARC','NUMPOINT',1,nPoint)
      call getrec(1,'JOBARC','NATOMS',1,nAtom)
      nSize = 3*nAtom

c   o set action flags
      bUpdate = (args(1)(1:3).eq.'upd')
      bDump   = (args(1)(1:4).eq.'dump').or.(args(1)(1:6).eq.'updump')
      bAnGrad = iflags2(138).eq.1

c   o initialize the icore free space pointer
      lFree = i0

      if (bUpdate.or.bDump) then

c      o find last displacement which was calculated
         call getrec(1,'JOBARC','FDCALCTP',nPoint,icore(lFree))
         iOff = nPoint-1
         do while ((icore(lFree+iOff).ge.0).and.(iOff.ge.0))
            iOff = iOff - 1
         end do
         iLast = 1+iOff
         if (iLast.eq.0) return

c      o print the reference energy
         if (bDump) then
            call getrec(1,'JOBARC','REFENERG',iintfp,dTmp)
            if (dTmp.ne.0.d0) call ja_dump('REFENERG',iintfp,1,dTmp)
         end if

c      o print the energy vector
         call getrec(1,'JOBARC','ENGPOINT',iintfp*nPoint,icore(lFree))
         if (bUpdate) then
            if (iflags(87).eq.0) then
               call getrec(1,'JOBARC','TOTENERG',iintfp,
     &                     icore(lFree+iintfp*iOff))
            else
               call getrec(1,'JOBARC','TOTENER2',iintfp,
     &                     icore(lFree+iintfp*iOff))
            end if
            call putrec(1,'JOBARC','ENGPOINT',iintfp*iLast,icore(lFree))
         end if
         if (bDump) then
            call ja_dump('ENGPOINT',iintfp,nPoint,icore(lFree))
         end if

         if (bAnGrad) then

            if (bUpdate) then
c            o determine the scope of symmetry
               call getrec(-1,'JOBARC','DANGERUS',1,iTmp)
               if (iflags(79).eq.0.and.iTmp.eq.0) then
                  szType='FULL'
               else
                  szType='COMP'
               end if
c            o assign indices
               lGrdInt = lFree
               lFree = lFree + iintfp*nSize
               lGrdXYZ = lFree
               lFree = lFree + iintfp*nSize
               lGeom = lFree
               lFree = lFree + iintfp*nSize
               lIMAP = lFree
               lFree = lFree + nAtom + iand(nAtom,1)
               iTmp = (icrsiz+i0-lFree)/iintfp
c            o load the original coordinates for the last displacement
               call getrec(1,'JOBARC','FDCOORDS',iintfp*nSize*iLast,
     &                     icore(lGeom))
               call c_memmove(icore(lGeom),
     &                        icore(lGeom+iintfp*nSize*iOff),
     &                        ifltln*nSize)
c            o load and rotate the derivatives
               call getgrd(nAtom,szType,icore(lGeom),
     &                     icore(lGrdXYZ),icore(lGrdInt),
     &                     dDipXYZ,dPolXYZ,icore(lIMAP),
     &                     icore(lFree),iTmp,.false.)
c            o reclaim scratch after GRDINT
               lFree = lGrdXYZ
            end if

c         o energy gradient
            call getrec(1,'JOBARC',
     &                  'GRDPOINT',iintfp*nSize*nPoint,icore(lFree))
            if (bUpdate) then
               call dcopy(nSize,icore(lGrdInt),1,
     &                          icore(lFree+iintfp*nSize*iOff),1)
               call putrec(1,'JOBARC',
     &                     'GRDPOINT',iintfp*nSize*iLast,icore(lFree))
            end if
            if (bDump) then
               Call ja_dump('GRDPOINT',iintfp*3,nAtom*nPoint,
     &                      icore(lFree))
            end if

c         o dipole
            call getrec(1,'JOBARC',
     &                  'DIPPOINT',iintfp*3*nPoint,icore(lFree))
            if (bUpdate) then
               call dcopy(3,dDipXYZ,1,icore(lFree+iintfp*3*iOff),1)
               call putrec(1,'JOBARC',
     &                     'DIPPOINT',iintfp*3*iLast,icore(lFree))
            end if
            if (bDump) then
               call ja_dump('DIPPOINT',iintfp*3,nPoint,
     &                      icore(lFree))
            end if

c         o polarizability
            call getrec(1,'JOBARC',
     &                  'POLPOINT',iintfp*9*nPoint,icore(lFree))
            if (bUpdate) then
               call dcopy(9,dPolXYZ,1,icore(lFree+iintfp*9*iOff),1)
               call putrec(1,'JOBARC',
     &                     'POLPOINT',iintfp*9*iLast,icore(lFree))
            end if
            if (bDump) then
               call ja_dump('POLPOINT',iintfp*3,3*nPoint,
     &                      icore(lFree))
            end if

c        end if (bAnGrad)
         end if

      else if (args(1)(1:4).eq.'load') then

c      o open the record file
         iLast = linblnk(args(2))
         inquire(file=args(2)(1:iLast),exist=bExist)
         if (bExist) then
            open(unit=10,file=args(2)(1:iLast),form='FORMATTED')
            rewind(10)
         else
            print *, '@PARFD: The record file does not exist.'
            print *, '        file = "',args(2)(1:iLast),'"'
            return
         end if

         bNotDone = .true.
         do while (bNotDone)

c         o read the record name
            read(10,fmt='(a)',end=100) szRecName

c         o branch on the record
            if      (szRecName(1:8).eq.'REFENERG') then
               call ja_load_dbl(10,szRecName(1:8),
     &                          1,1,icore(lFree))
            else if (szRecName(1:8).eq.'ENGPOINT') then
               call ja_load_dbl(10,szRecName(1:8),
     &                          1,nPoint,icore(lFree))
            else if (szRecName(1:8).eq.'GRDPOINT') then
               call ja_load_dbl(10,szRecName(1:8),
     &                          3,nAtom*nPoint,icore(lFree))
            else if (szRecName(1:8).eq.'DIPPOINT') then
               call ja_load_dbl(10,szRecName(1:8),
     &                          3,nPoint,icore(lFree))
            else if (szRecName(1:8).eq.'POLPOINT') then
               call ja_load_dbl(10,szRecName(1:8),
     &                          3,3*nPoint,icore(lFree))
            else
               bNotDone = .false.
            end if

         end do
 100     continue

c      o close the record file
         close(unit=10,status='KEEP')

      else

         print *, '@PARFD: The parfd module requires 1 or 2 arguments'
         print '()'
         print *, '  xa2proc parfd update'
         print *, '  xa2proc parfd (up)dump'
         print *, '  xa2proc parfd load <file>'
         call aces_exit(1)

      end if

      return
c     end subroutine parfd
      end

c ----------------------------------------------------------------------

      subroutine ja_dump(szRecName,nRows,nCols,icore)
      implicit none

      character*8 szRecName
      integer nRows, nCols
      integer icore(nRows,nCols)

      integer iRow, iCol

      if (nRows.eq.0.or.nCols.eq.0) return

      print '(a)', szRecName
      do iCol = 1, nCols
         print *, (icore(iRow,iCol),iRow=1,nRows)
      end do

      return
      end

c ----------------------------------------------------------------------

      subroutine ja_load_dbl(iUnit,szRecName,nRows,nCols,icore)
      implicit none

      integer iUnit
      character*8 szRecName
      integer nRows, nCols
      integer icore(*)

      integer iOff, iRow, iCol



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end




      iOff = 0
      do iCol = 1, nCols
         read(10,fmt=*) (icore(iRow+iOff),iRow=1,iintfp*nRows)
         iOff = iOff + iintfp*nRows
      end do
C      print *, 'DEBUG: icore = ',(icore(iOff),iOff=1,iintfp*nRows*nCols)

      call getrec(1,'JOBARC',szRecName,iintfp*nRows*nCols,icore(1+iOff))
      call parfd_repl_nz(nRows*nCols,icore(1),icore(1+iOff))
      call putrec(1,'JOBARC',szRecName,iintfp*nRows*nCols,icore(1+iOff))

      return
      end

c ----------------------------------------------------------------------

      subroutine parfd_repl_nz(iNum,dSrc,dDest)
      implicit none
      integer iNum
      double precision dSrc(*), dDest(*)
      integer i
      do i = 1, iNum
         if (dSrc(i).ne.0.d0) dDest(i) = dSrc(i)
      end do
      return
      end

