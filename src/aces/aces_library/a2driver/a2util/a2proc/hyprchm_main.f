
c This program writes a drives creation of HYPERCHEM  format file containing 
c the coordinates, frequencies, normal mode vectors and intensities.
c
c for frequencies NORMCO
c 
c Using Ken Wilson's March 1998, Marshall Cory & Ajith Perera, 10/04.

      subroutine hyprchm_main
      implicit double precision (a-h,o-z)
      logical tbohr,hfcrap,linear
C


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end














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



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end













































































































































































































































































































































































































































































C
      maxcor=icrsiz
C
      call getrec(20,'JOBARC','NATOMS',1,natoms)

      iatchrg  = i0
      icoord   = iatchrg  + natoms+mod(natoms,2)
      ifuchrg  = icoord   + 3*natoms*iintfp
      ifucoord = ifuchrg  + natoms+mod(natoms,2)
      ifreqco  = ifucoord + 3*natoms*iintfp
      ifreq    = ifreqco  + 3*natoms*iintfp
      inormmd  = ifreq    + 3*natoms*iintfp
      intnsty  = inormmd  + 9*natoms*natoms*iintfp
      iord     = intnsty  + 3*natoms*iintfp
      itag     = iord     + 3*natoms
      inext    = itag     + 10*natoms
C
      if(inext-i0.gt.maxcor)call insmem('hyprchm-main',inext-i0,
     &                                   maxcor)

      call hyprchm_rd_geomvib(natoms,icore(iatchrg),icore(icoord),
     &                        icore(ifuchrg),icore(ifucoord),
     &                        icore(ifreqco),icore(ifreq),
     &                        icore(inormmd),icore(intnsty),
     &                        icore(iord),hfcrap,linear)


      tbohr = .true.
      if (iflags(54) .ne.0) hfcrap = .true.

      if(natoms .gt. 1) call wrt_hyprchm( natoms, icore(icoord), 
     .                  icore(ifreq), icore(inormmd), icore(intnsty),
     .                  icore(iatchrg),tbohr,hfcrap,linear,icore(iord),
     .                  icore(itag))

      print *, '@HYPERCHEM: successfully created interface files'

      return
      end

