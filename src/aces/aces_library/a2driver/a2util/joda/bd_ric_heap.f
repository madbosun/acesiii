
c This routine defines the default values for the common block in
c include/ric_heap.com.

      blockdata bd_ric_heap
      implicit none
c ric_heap.com : begin

c This common block contains the heap address and array indices for
c processing RICs. New additions must be initialized to 1 in
c bd_ric_heap.F and set in init_ric_heap.F.


      double precision dRICHeap(1)
      integer          z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv

      common /ric_heap_com/ dRICHeap,
     &                      z_RICHeap, z_DerBMat, z_BMat, z_GMat,
     &                      z_BTGInv
      save   /ric_heap_com/

c ric_heap.com : end
      data dRICHeap /0.d0/
      data z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv /5*1/
      end

