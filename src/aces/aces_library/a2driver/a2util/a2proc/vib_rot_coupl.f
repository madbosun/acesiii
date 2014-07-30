













































































































































































































































































































































































































































































































































      subroutine vib_rot_coupl
C
      implicit double precision (a-h,o-z)
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
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
C-----------------------------------------------------------------------
      if (iflags(11).eq.0) iuhf = 0
      maxcor=icrsiz
C
      Call Getrec(20,'JOBARC','NREALATM',1,Nreals)
      Call Getrec(20,'JOBARC','ZMATATMS',1,Natoms)
      Call Getrec(20, 'JOBARC', 'LINEAR', 1, Iamlinear)
C
      If (Iamlinear .EQ. 1) Then
         Nxm6 = 3*Nreals - 5
      Else
         Nxm6 = 3*Nreals - 6
      Endif
      Nx  = 3*Nreals
C
      iatmMass = i0
      icoord   = iatmMass + 3*Natoms*iintfp 
      iBmat    = icoord   + 3*NAtoms*iintfp
      iGmat    = iBmat    + Nx*NXm6*iintfp
      ihess    = iGmat    + Nxm6*Nxm6*iintfp
      ieig_vec = ihess    + Nx*Nx*iintfp
      ivib_vec = ieig_vec + Nx*Nx*iintfp
      ibtmp    = ivib_vec + Nx*Nxm6*iintfp
      itmp1    = ibtmp    + Nx*Nxm6*iintfp 
      itmp2    = itmp1    + Nx*Nx*iintfp 
      immatx   = Itmp2    + Nx*Nx*iintfp
      immaty   = Immatx   + Nx*Nx*iintfp
      immatz   = Immaty   + Nx*Nx*iintfp
      imap     = Immatz   + Nx*Nx*iintfp
      inext    = Imap     + Natoms
C
      if(inext-i0.gt.maxcor)call insmem('vib_rot_coupl',inext-i0,
     &                      maxcor)

      call coriolis_cc(natoms,nreals,Nx,Nxm6,icore(iatmMass),
     &                 icore(icoord),icore(ibmat),icore(igmat),
     &                 icore(ihess),icore(ieig_vec),icore(ivib_vec),
     &                 icore(ibtmp),icore(itmp1),icore(itmp2),
     &                 icore(immatx),icore(immaty),icore(immatz),
     &                 icore(imap),iuhf)
c ----------------------------------------------------------------------

      return
      end

