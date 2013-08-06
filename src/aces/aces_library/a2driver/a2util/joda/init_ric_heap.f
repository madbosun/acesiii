
c This routine allocates a heap for processing RICs and creates pointers
c into the heap for various arrays.

      subroutine init_ric_heap
      implicit none



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






c ric_heap.com : begin

c This common block contains the heap address and array indices for
c processing RICs. New additions must be initialized to 1 in
c bd_ric_heap.F and set in init_ric_heap.F.

      external bd_ric_heap

      double precision dRICHeap(1)
      integer          z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv

      common /ric_heap_com/ dRICHeap,
     &                      z_RICHeap, z_DerBMat, z_BMat, z_GMat,
     &                      z_BTGInv
      save   /ric_heap_com/

c ric_heap.com : end
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

      integer iNeed, ndx
      integer s_DerBMat, s_BMat, s_GMat, s_BTGInv

      integer vp(2)
      data vp /0,0/

c   o don't allocate memory more than once
      call c_memmove(vp,dRICHeap,ifltln)
      if (vp(1).ne.0.or.vp(2).ne.0) return

      iNeed = 0

c   o DerBMat: The derivative of the B-matrix for transforming the
c              Hessian between RICs and Cartesians.
      s_DerBMat = 9*mxatms*mxatms*maxredunco
      iNeed = iNeed + s_DerBMat

c   o BMat: The B-matrix for RIC/Cartesian transformations.
      s_BMat = 9*mxatms*mxatms
      iNeed = iNeed + s_BMat

c   o GMat: The G-matrix for RIC/Cartesian gradient transformations.
      s_GMat = 9*mxatms*mxatms
      iNeed = iNeed + s_GMat

c   o BTGInv: inv(trans(B)*G) required for transforming the Hessian.
      s_BTGInv = 3*mxatms*maxredunco
      iNeed = iNeed + s_BTGInv

      ndx = iintfp*iNeed
      call aces_malloc(ndx,dRICHeap,z_RICHeap)
      call c_memmove(vp,dRICHeap,ifltln)
      if (vp(1).eq.0.and.vp(2).eq.0) then
         print *, '@INIT_RIC_HEAP: Failed to allocate memory.'
         print *, '                need ',iNeed/1024/1024,' MB'
         call errex
      end if
      z_RICHeap = (z_RICHeap+iintfp-1)/iintfp
      ndx = z_RICHeap

c   o DerBMat location
      z_DerBMat = ndx
      ndx = ndx + s_DerBMat

c   o BMat location
      z_BMat = ndx
      ndx = ndx + s_BMat

c   o GMat location
      z_GMat = ndx
      ndx = ndx + s_GMat

c   o BTGInv location
      z_BTGInv = ndx
      ndx = ndx + s_BTGInv

c      print *, '@INIT_RIC_HEAP: RIC HEAP INDICES'
c      print *, '                z_RICHeap = ',z_RICHeap
c      print *, '                z_DerBMat = ',z_DerBMat
c      print *, '                z_BMat    = ',z_BMat
c      print *, '                z_GMat    = ',z_GMat
c      print *, '                z_BTGInv  = ',z_BTGInv

      return
c     end subroutine init_ric_heap
      end

